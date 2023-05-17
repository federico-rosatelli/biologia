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
	typeSearch := r.URL.Query().Get("type")
	Println(search)
	p, err := rt.db.FindTaxon(search, typeSearch)
	if err != nil {
		http.Error(w, err.Error(), http.StatusBadRequest)
		return
	}

	if errJson := json.NewEncoder(w).Encode(p); errJson != nil {
		http.Error(w, errJson.Error(), http.StatusBadRequest)
		return
	}
}
