package main

import (
	"errors"
	"flag"
	"fmt"
	"io"
	"os"
	"time"

	"github.com/ardanlabs/conf"
	"gopkg.in/yaml.v3"
)

type WebConfig struct {
	Config struct {
		Path      string `conf:"default:./conf/config.yml"`
		PathDebug string `conf:"default:./conf/configDebug.yml"`
	}
	WebAPI struct {
		ServerName    string `conf:"default:0.0.0.0"`
		ServerVersion string `conf:"default:a-0.0.1"`
		ServerPort    int    `conf:"default:3000"`
		ServerAddr    string `conf:"default:0.0.0.0:8080"`
		ServerConfig  struct {
			ReadTimeout     time.Duration `conf:"default:5s"`
			WriteTimeout    time.Duration `conf:"default:5s"`
			ShutdownTimeout time.Duration `conf:"default:5s"`
		}
	}
	Debug  bool
	DevRun bool
}

// loadConfig creates a WebConfig from flags and configuration file.
func loadConfig() (WebConfig, error) {
	var cfg WebConfig

	if err := conf.Parse(os.Args[1:], "CFG", &cfg); err != nil {
		if errors.Is(err, conf.ErrHelpWanted) {
			usage, err := conf.Usage("CFG", &cfg)
			if err != nil {
				return cfg, fmt.Errorf("Error generating config usage: %w", err)
			}
			fmt.Println(usage) //nolint:forbidigo
			return cfg, conf.ErrHelpWanted
		}
		return cfg, fmt.Errorf("parsing config: %w", err)
	}

	fp, err := os.Open(cfg.Config.Path)
	if err != nil && !os.IsNotExist(err) {
		return cfg, fmt.Errorf("Can't read the config file: %w", err)
	} else if err == nil {
		yamlFile, err := io.ReadAll(fp)
		if err != nil {
			return cfg, fmt.Errorf("can't read config file: %w", err)
		}
		err = yaml.Unmarshal(yamlFile, &cfg)
		if err != nil {
			return cfg, fmt.Errorf("can't unmarshal config file: %w", err)
		}
		_ = fp.Close()
	}
	boolFlag := flag.Bool("dev", false, "Database Usage for dev")
	flag.Parse()
	cfg.DevRun = !*boolFlag

	return cfg, nil
}
