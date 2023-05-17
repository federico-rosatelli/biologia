package api

import (
	"net/http"
)

func (rt *_router) Handler() http.Handler {
	// Register routes
	rt.router.GET("/", rt.Welcome)
	rt.router.GET("/taxonomy", rt.Taxonomy)
	rt.router.GET("/taxonomy1", rt.Taxon_Term)
	return rt.router
}
