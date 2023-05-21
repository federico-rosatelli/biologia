package database

import (
	errorM "GOServer/service/ErrManager"
	str "GOServer/service/structures"
	"context"

	"go.mongodb.org/mongo-driver/bson"
	"go.mongodb.org/mongo-driver/mongo/options"
)

func (db *appDB) FindProteinsId(taxonId string) (str.TableBasic, errorM.Errors) {
	var nBasic str.TableBasic
	proj := options.FindOne().SetProjection(bson.D{{Key: "Nucleotides", Value: 0}})
	errM := db.table_basic.FindOne(context.TODO(), bson.D{{Key: "TaxId", Value: taxonId}}, proj).Decode(&nBasic)
	if errM != nil {
		return nBasic, errorM.NewError(errM.Error(), errorM.StatusBadRequest)
	}
	return nBasic, nil
}

func (db *appDB) FindProteinTableByLocus(taxId string, locus string) (str.TableComplete, errorM.Errors) {
	var protein str.TableComplete
	proj := options.FindOne().SetProjection(bson.D{{Key: "Proteins", Value: bson.M{"$elemMatch": bson.M{"GBSeq_locus": locus}}}})
	errM := db.table_complete.FindOne(context.TODO(), bson.M{"TaxId": taxId, "Proteins.GBSeq_locus": locus}, proj).Decode(&protein)
	if errM != nil {
		return protein, errorM.NewError(errM.Error(), errorM.StatusBadRequest)
	}
	return protein, nil

}

func (db *appDB) FindProteinByLocus(locus string) (str.Protein, errorM.Errors) {
	var protein str.Protein
	errM := db.protein_data.FindOne(context.TODO(), bson.D{{Key: "GBSeq_locus", Value: locus}}).Decode(&protein)
	if errM != nil {
		return protein, errorM.NewError(errM.Error(), errorM.StatusBadRequest)
	}
	return protein, nil

}
