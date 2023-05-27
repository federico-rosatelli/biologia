#!/usr/bin/env sh

docker build -t bio-frontend:latest -f Dockerfile.frontend .
docker build -t bio-backend:latest -f Dockerfile.backend .