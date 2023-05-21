package api

import (
	"GOServer/service/ErrManager"
	str "GOServer/service/structures"
)

func (rt *_router) GetProteins(taxId string) (str.TableBasic, ErrManager.Errors) {
	nBasic, err := rt.db.FindProteinsId(taxId)
	if err != nil {
		return nBasic, err
	}
	nu := myNucleotideBasic(nBasic)

	return nBasic, nu.auxProteinBasic()
}

func (nucl myNucleotideBasic) auxProteinBasic() ErrManager.Errors {
	if len(nucl.Proteins) == 0 {
		return ErrManager.NewError("No Content To Show", ErrManager.StatusServiceUnavailable)
	}
	if nucl.ScientificName == "" {
		return ErrManager.NewError("Content without Name", ErrManager.StatusBadRequest)
	}
	return nil
}

func (rt *_router) GetProteinFromTable(taxId string, locus string) (str.TableComplete, ErrManager.Errors) {
	var tableProtein str.TableComplete
	tableProtein, err := rt.db.FindProteinTableByLocus(taxId, locus)
	if err != nil {
		protein, errM := rt.db.FindProteinByLocus(locus)
		if errM != nil {
			return tableProtein, errM
		}
		tableProtein.Proteins = append(tableProtein.Proteins, protein)
	}
	return tableProtein, nil
}
