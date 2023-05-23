package api

import (
	"encoding/json"
	"net/http"

	"github.com/julienschmidt/httprouter"
)

func (rt *_router) Welcome(w http.ResponseWriter, r *http.Request, ps httprouter.Params) {
	w.Header().Set("content-type", "application/json")
	if errJson := json.NewEncoder(w).Encode("OK"); errJson != nil {
		http.Error(w, errJson.Error(), http.StatusBadRequest)
		return
	}
}

// Da rinonimare
func (rt *_router) Taxonomy(w http.ResponseWriter, r *http.Request, ps httprouter.Params) {
	w.Header().Set("content-type", "application/json")

	search := r.URL.Query().Get("search")
	Println(search)
	p, err := rt.db.FindTaxon(search)
	if err != nil {
		http.Error(w, err.Error(), http.StatusBadRequest)
		return
	}

	if errJson := json.NewEncoder(w).Encode(p); errJson != nil {
		http.Error(w, errJson.Error(), http.StatusBadRequest)
		return
	}
}

func (rt *_router) TaxonomyTree(w http.ResponseWriter, r *http.Request, ps httprouter.Params) {
	w.Header().Set("content-type", "application/json")

	search := r.URL.Query().Get("search")
	Println(search)
	p, err := rt.db.FindTaxonTree(search)
	if err != nil {
		http.Error(w, err.Error(), http.StatusBadRequest)
		return
	}

	if errJson := json.NewEncoder(w).Encode(p); errJson != nil {
		http.Error(w, errJson.Error(), http.StatusBadRequest)
		return
	}
}

func (rt *_router) Taxon_Term(w http.ResponseWriter, r *http.Request, ps httprouter.Params) {
	w.Header().Set("content-type", "application/json")

	search := r.URL.Query().Get("search")
	typeSearch := r.URL.Query().Get("type")
	location := r.URL.Query().Get("location")
	product := r.URL.Query().Get("product")
	p, err := rt.GetOrganisms(search, typeSearch, location, product)
	if err != nil {
		http.Error(w, err.Error(), err.Type())
		return
	}

	if errJson := json.NewEncoder(w).Encode(p); errJson != nil {
		http.Error(w, errJson.Error(), http.StatusBadRequest)
		return
	}
}

func (rt *_router) Get_Nucleotides_Id(w http.ResponseWriter, r *http.Request, ps httprouter.Params) {
	w.Header().Set("content-type", "application/json")

	id := ps.ByName("id")
	p, err := rt.GetNucleotides(id)
	if err != nil {
		http.Error(w, err.Error(), err.Type())
		return
	}

	if errJson := json.NewEncoder(w).Encode(p); errJson != nil {
		http.Error(w, errJson.Error(), http.StatusBadRequest)
		return
	}
}

func (rt *_router) Get_Nucleotide_From_Locus(w http.ResponseWriter, r *http.Request, ps httprouter.Params) {
	w.Header().Set("content-type", "application/json")

	id := ps.ByName(("id"))
	locus := ps.ByName("locus")
	table, err := rt.GetNucleotideFromTable(id, locus)
	if err != nil {
		http.Error(w, err.Error(), err.Type())
		return
	}

	if errJson := json.NewEncoder(w).Encode(table); errJson != nil {
		http.Error(w, errJson.Error(), http.StatusBadRequest)
		return
	}

}

func (rt *_router) Get_Proteins_Id(w http.ResponseWriter, r *http.Request, ps httprouter.Params) {
	w.Header().Set("content-type", "application/json")

	id := ps.ByName("id")
	p, err := rt.GetProteins(id)
	if err != nil {
		http.Error(w, err.Error(), err.Type())
		return
	}

	if errJson := json.NewEncoder(w).Encode(p); errJson != nil {
		http.Error(w, errJson.Error(), http.StatusBadRequest)
		return
	}
}

func (rt *_router) Get_Protein_From_Locus(w http.ResponseWriter, r *http.Request, ps httprouter.Params) {
	w.Header().Set("content-type", "application/json")

	id := ps.ByName(("id"))
	locus := ps.ByName("locus")
	table, err := rt.GetProteinFromTable(id, locus)
	if err != nil {
		http.Error(w, err.Error(), err.Type())
		return
	}

	if errJson := json.NewEncoder(w).Encode(table); errJson != nil {
		http.Error(w, errJson.Error(), http.StatusBadRequest)
		return
	}

}
