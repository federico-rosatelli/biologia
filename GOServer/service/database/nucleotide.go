package database

import (
	errorM "GOServer/service/ErrManager"
	str "GOServer/service/structures"
	"context"

	"go.mongodb.org/mongo-driver/bson"
	"go.mongodb.org/mongo-driver/mongo/options"
)

func (db *appDB) FindNucleotidesId(taxonId string) (str.TableBasic, errorM.Errors) {
	var nBasic str.TableBasic
	var nCompl str.TableBasic
	proj := options.FindOne().SetProjection(bson.D{{Key: "Proteins", Value: 0}})
	errM := db.table_basic.FindOne(context.TODO(), bson.D{{Key: "TaxId", Value: taxonId}}, proj).Decode(&nBasic)
	if errM != nil {
		return nBasic, errorM.NewError(errM.Error(), errorM.StatusBadRequest)
	}
	proj = options.FindOne().SetProjection(bson.D{{Key: "Products", Value: 1}})
	errM = db.table_complete.FindOne(context.TODO(), bson.D{{Key: "TaxId", Value: taxonId}}, proj).Decode(&nCompl)
	if errM != nil {
		return nBasic, errorM.NewError(errM.Error(), errorM.StatusBadRequest)
	}
	nBasic.Products = nCompl.Products
	return nBasic, nil
}

func (db *appDB) FindNucleotideTableByLocus(taxId string, locus string) (str.TableComplete, errorM.Errors) {
	var nucleo str.TableComplete
	proj := options.FindOne().SetProjection(bson.D{{Key: "Nucleotides", Value: bson.M{"$elemMatch": bson.M{"GBSeq_locus": locus}}}})
	errM := db.table_complete.FindOne(context.TODO(), bson.M{"TaxId": taxId, "Nucleotides.GBSeq_locus": locus}, proj).Decode(&nucleo)
	if errM != nil {
		return nucleo, errorM.NewError(errM.Error(), errorM.StatusBadRequest)
	}
	return nucleo, nil

}

func (db *appDB) FindNucleotideByLocus(locus string) (str.Nucleotide, errorM.Errors) {
	var nucleo str.Nucleotide
	errM := db.nucleotide_data.FindOne(context.TODO(), bson.D{{Key: "GBSeq_locus", Value: locus}}).Decode(&nucleo)
	if errM != nil {
		return nucleo, errorM.NewError(errM.Error(), errorM.StatusBadRequest)
	}
	return nucleo, nil

}
