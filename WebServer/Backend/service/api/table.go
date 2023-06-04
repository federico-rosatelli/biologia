package api

import (
	"GOServer/service/ErrManager"
	str "GOServer/service/structures"
	"strings"
)

type myOrganismTable []str.OrganismTable

func (rt *_router) GetOrganisms(search string, nType string, location string, product string) ([]str.OrganismTable, ErrManager.Errors) {
	var p []str.OrganismTable
	if nType == "" {
		return p, ErrManager.NewError("Search & Type Cannot Be Empty", ErrManager.StatusBadRequest)
	}
	if (search == "" || strings.TrimSpace(search) == "") && (location == "" || strings.TrimSpace(location) == "") && (product == "" || strings.TrimSpace(product) == "") {
		return p, ErrManager.NewError("Error Empty String", ErrManager.StatusBadRequest)
	}
	var listTaxId []string
	if product != "" {
		listTaxId_aux, err := rt.db.TableOrganismByProduct(product)
		if err != nil {
			return p, err
		}
		if len(listTaxId_aux) == 0 {
			return p, ErrManager.NewError("Product Not Found", ErrManager.StatusBadRequest)
		}
		listTaxId = append(listTaxId, listTaxId_aux...)
	}
	if location != "" {
		listTaxId_aux, err := rt.db.TableOrganismByLocation(location)
		if err != nil {
			return p, err
		}
		if len(listTaxId_aux) == 0 {
			return p, ErrManager.NewError("Location Not Found", ErrManager.StatusBadRequest)
		}
		listTaxId = append(listTaxId, listTaxId_aux...)
	}

	p, err := rt.db.TableOrganism(search, nType, listTaxId)
	if err != nil {
		return p, err
	}
	table := myOrganismTable(p)
	return p, table.auxOrganismTable()
}

func (table myOrganismTable) auxOrganismTable() ErrManager.Errors {
	if len(table) == 0 {
		return ErrManager.NewError("No Content Found", ErrManager.StatusBadRequest)
	}
	for _, t := range table {
		if t.ScientificName == "" || t.TaxId == "" {
			return ErrManager.NewError("Contentent Must have ScientificName&TaxId. "+t.ScientificName, ErrManager.StatusBadRequest)
		}
	}
	return nil
}
