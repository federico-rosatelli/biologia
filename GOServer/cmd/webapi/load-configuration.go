package main

import (
	"os"
	"time"

	"gopkg.in/yaml.v3"
)

type WebServerConfig struct {
	Name    string `yaml:"Name"`
	Version string `yaml:"Version"`
	Server  struct {
		ServerPort    int    `yaml:"ServerPort"`
		ServerIp      string `yaml:"ServerIp"`
		ServerAddress string `yaml:"ServerAddress"`
		ServerConfig  struct {
			ReadTimeout     time.Duration `yaml:"ReadTimeout"`
			WriteTimeout    time.Duration `yaml:"WriteTimeout"`
			ShutdownTimeout time.Duration `yaml:"ShutdownTimeout"`
		}
	} `yaml:"Server"`
}

// loadConfig creates a WebConfig from flags and configuration file.
func loadConfig() (WebServerConfig, error) {
	var wsc WebServerConfig

	file, err := os.Open("./conf/config.yml")
	if err != nil {
		return wsc, err
	}
	defer file.Close()

	d := yaml.NewDecoder(file)
	if err := d.Decode(&wsc); err != nil {
		return wsc, err
	}

	return wsc, nil
}
