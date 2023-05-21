package api

import (
	"net/http"
)

func (rt *_router) Handler() http.Handler {
	// Register routes
	rt.router.GET("/", rt.Welcome)
	rt.router.GET("/taxonomy", rt.Taxonomy)
	rt.router.GET("/taxon_term", rt.Taxon_Term)
	rt.router.GET("/organism/:id/nucleotides", rt.Get_Nucleotides_Id)
	rt.router.GET("/organism/:id/nucleotides/:locus", rt.Get_Nucleotide_From_Locus)
	rt.router.GET("/organism/:id/proteins", rt.Get_Proteins_Id)
	rt.router.GET("/organism/:id/proteins/:locus", rt.Get_Protein_From_Locus)
	return rt.router
}
