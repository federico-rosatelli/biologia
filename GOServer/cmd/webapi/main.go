package main

import (
	"context"
	"errors"
	"fmt"
	"net/http"
	"os"
	"os/signal"
	"syscall"

	"GOServer/service/api"
	"GOServer/service/database"

	"github.com/ardanlabs/conf"
	"github.com/sirupsen/logrus"
	"go.mongodb.org/mongo-driver/mongo"
	"go.mongodb.org/mongo-driver/mongo/options"
)

func main() {
	if err := run(); err != nil {
		_, _ = fmt.Fprintln(os.Stderr, "Error running the app: ", err)
	}
}

func run() error {
	cfg, err := loadConfig()
	if err != nil {
		if errors.Is(err, conf.ErrHelpWanted) {
			return nil
		}
		return err
	}

	logger := logrus.New()
	logger.Infof("Project \t%s", cfg.Name)
	logger.Infof("Version \t%s", cfg.Version)
	logger.Infof("Server:%s Starting on Port: %d", cfg.Server.ServerIp, cfg.Server.ServerPort)
	var db database.AppDatabase

	logger.Println("Initializing Database Support")
	dbconn := options.Client().ApplyURI("mongodb://localhost:27017/")
	client, err := mongo.Connect(context.TODO(), dbconn)
	if err != nil {
		logger.WithError(err).Error("error connetting to mongo DB")
		return fmt.Errorf("opening mongoDb: %w", err)
	}
	defer func() {
		logger.Debug("database stopping")
		_ = client.Disconnect(context.TODO())
	}()
	db, err = database.InitDatabase(client)
	if err != nil {
		logger.WithError(err).Error("error creating mongo DB")
		return fmt.Errorf("opening mongoDb: %w", err)
	}

	shutdown := make(chan os.Signal, 1)
	signal.Notify(shutdown, os.Interrupt, syscall.SIGTERM)

	serverErrors := make(chan error, 1)

	apirouter, err := api.NewServer(api.Config{
		Logger:   logger,
		Database: db,
	})

	if err != nil {
		logger.WithError(err).Error("error creating the API server instance")
		return fmt.Errorf("creating the API server instance: %w", err)
	}
	router := apirouter.Handler()

	router, err = registerWebUI(router)
	if err != nil {
		logger.WithError(err).Error("error registering web UI handler")
		return fmt.Errorf("registering web UI handler: %w", err)
	}

	router = applyCORSHandler(router)

	apiserver := http.Server{
		Addr:              cfg.Server.ServerAddress,
		Handler:           router,
		ReadTimeout:       cfg.Server.ServerConfig.ReadTimeout,
		ReadHeaderTimeout: cfg.Server.ServerConfig.ReadTimeout,
		WriteTimeout:      cfg.Server.ServerConfig.WriteTimeout,
	}

	go func() {
		logger.Infof("API listening on port: %d", cfg.Server.ServerPort)
		serverErrors <- apiserver.ListenAndServe()
		logger.Infof("stopping API server")
	}()

	select {
	case err := <-serverErrors:

		return fmt.Errorf("server error: %w", err)

	case sig := <-shutdown:
		logger.Infof("signal %v received, start shutdown", sig)

		err := apirouter.Close()
		if err != nil {
			logger.WithError(err).Warning("graceful shutdown of apirouter error")
		}

		ctx, cancel := context.WithTimeout(context.Background(), cfg.Server.ServerConfig.ShutdownTimeout)
		defer cancel()

		err = apiserver.Shutdown(ctx)
		if err != nil {
			logger.WithError(err).Warning("error during graceful shutdown of HTTP server")
			err = apiserver.Close()
		}

		switch {
		case sig == syscall.SIGSTOP:
			return errors.New("integrity issue caused shutdown")
		case err != nil:
			return fmt.Errorf("could not stop server gracefully: %w", err)
		}
	}

	return nil
}
