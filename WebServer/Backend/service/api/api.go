package api

import (
	"GOServer/service/database"
	"errors"
	"net/http"
	"time"

	"github.com/julienschmidt/httprouter"
	"github.com/sirupsen/logrus"
)

type Router interface {
	Handler() http.Handler

	Close() error
}

type Config struct {
	Logger   logrus.FieldLogger
	Database database.AppDatabase
}

type _router struct {
	router *httprouter.Router

	baseLogger logrus.FieldLogger

	db database.AppDatabase
}

func NewServer(cfg Config) (Router, error) {
	if cfg.Logger == nil {
		return nil, errors.New("logger is required")
	}
	router := httprouter.New()
	router.RedirectTrailingSlash = false
	router.RedirectFixedPath = false

	// if database

	return &_router{
		router:     router,
		baseLogger: cfg.Logger,
		db:         cfg.Database,
	}, nil
}

func (rt *_router) Close() error {
	return nil
}

func Println(log any, params ...int) {
	logger := logrus.New()
	level := 4
	if len(params) != 0 {
		level = params[0]
	}
	logger.SetLevel(logrus.Level(level))
	if level == 3 {
		logger.Log(logrus.Level(level), time.Now().String())
	}
	logger.Log(logrus.Level(level), log)
}
