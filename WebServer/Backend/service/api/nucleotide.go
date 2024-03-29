package api

import (
	"GOServer/service/ErrManager"
	str "GOServer/service/structures"
)

type myTableBasic str.TableBasic

func (rt *_router) GetNucleotides(taxId string) (str.TableBasic, ErrManager.Errors) {
	nBasic, err := rt.db.FindNucleotidesId(taxId)
	if err != nil {
		return nBasic, err
	}
	nu := myTableBasic(nBasic)

	return nBasic, nu.auxNucleoBasic()
}

func (nucl myTableBasic) auxNucleoBasic() ErrManager.Errors {
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
