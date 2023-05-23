package api

import (
	"GOServer/service/ErrManager"
	str "GOServer/service/structures"
	"strings"
)

func (rt *_router) GetTaxonomy(taxId string) (str.Taxonomy, ErrManager.Errors) {
	var taxonomy str.Taxonomy
	taxonomy, err := rt.db.FindTaxon(taxId)
	if err != nil {
		return taxonomy, err
	}
	return taxonomy, nil
}

func (rt *_router) GetTaxonomyTree(search string, nType string) (str.TaxonomyTree, ErrManager.Errors) {
	var tree str.TaxonomyTree
	if nType == "" {
		return tree, ErrManager.NewError("Search & Type Cannot Be Empty", ErrManager.StatusBadRequest)
	}
	if search == "" || strings.TrimSpace(search) == "" {
		return tree, ErrManager.NewError("Error Empty String", ErrManager.StatusBadRequest)
	}
	tree, err := rt.db.FindTaxonTree(search, nType)
	if err != nil {
		return tree, err
	}
	return tree, nil
}
