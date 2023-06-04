package main

import (
	"net/http"

	"github.com/gorilla/handlers"
)

// applyCORSHandler applies a CORS policy to the router. CORS stands for Cross-Origin Resource Sharing: it's a security
// feature present in web browsers that blocks JavaScript requests going across different domains if not specified in a
// policy. This function sends the policy of this API server.
func applyCORSHandler(h http.Handler) http.Handler {
	return handlers.CORS(
		handlers.AllowedMethods([]string{"GET", "POST", "OPTIONS", "DELETE", "PUT"}),
		handlers.AllowedOrigins([]string{"*"}),
		handlers.AllowedHeaders([]string{"Content-Type", "Token", "X-Requested-With", "Content-Type", "Origin", "Cache-Control", "Pragma", "Authorization", "Accept", "Accept-Encoding"}),
		handlers.ExposedHeaders([]string{"Token"}),
		handlers.MaxAge(10),
		handlers.AllowCredentials(),
	)(h)
}
