package api

import (
	"GOServer/service/ErrManager"
	str "GOServer/service/structures"
)

type myNucleotideBasic str.TableBasic

type myOrganismTable []str.OrganismTable

func (rt *_router) GetOrganisms(search string, nType string) ([]str.OrganismTable, ErrManager.Errors) {
	var p []str.OrganismTable
	if search == "" || search == " " || nType == "" {
		return p, ErrManager.NewError("Search & Type Cannot Be Empty", ErrManager.StatusBadRequest)
	}
	p, err := rt.db.TableOrganism(search, nType)
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

func (rt *_router) GetNucleotides(taxId string) (str.TableBasic, ErrManager.Errors) {
	nBasic, err := rt.db.FindNucleotidesId(taxId)
	if err != nil {
		return nBasic, err
	}
	nu := myNucleotideBasic(nBasic)

	return nBasic, nu.auxNucleoBasic()
}

func (nucl myNucleotideBasic) auxNucleoBasic() ErrManager.Errors {
	if len(nucl.Nucleotides) == 0 {
		return ErrManager.NewError("No Content To Show", ErrManager.StatusServiceUnavailable)
	}
	if nucl.ScientificName == "" {
		return ErrManager.NewError("Content without Name", ErrManager.StatusBadRequest)
	}
	return nil
}

func (rt *_router) GetNucleotideFromTable(taxId string, locus string) (str.TableComplete, ErrManager.Errors) {
	var tableNucleo str.TableComplete
	tableNucleo, err := rt.db.FindNucleotideTableByLocus(taxId, locus)
	if err != nil {
		nucleo, errM := rt.db.FindNucleotideByLocus(locus)
		if errM != nil {
			return tableNucleo, errM
		}
		tableNucleo.Nucleotides = append(tableNucleo.Nucleotides, nucleo)
	}
	return tableNucleo, nil
}
