package database

import (
	errorM "GOServer/service/ErrManager"
	str "GOServer/service/structures"
	"context"

	"go.mongodb.org/mongo-driver/bson"
	"go.mongodb.org/mongo-driver/mongo"
	"go.mongodb.org/mongo-driver/mongo/options"
)

func (db *appDB) TableOrganism(search string, typeS string) ([]str.OrganismTable, errorM.Errors) {
	var orgTable []str.OrganismTable
	filter := bson.M{}
	if typeS == "id" {
		filter["TaxId"] = search
	} else if typeS == "scientific_name" {
		filter["ScientificName"] = bson.D{{Key: "$regex", Value: search}}
	} else {
		return orgTable, errorM.NewError("Only id|scientific_name allowed", errorM.StatusBadRequest)
	}

	groupStage := bson.D{{Key: "$project", Value: bson.M{"QtyNucleotides": bson.M{"$size": "$Nucleotides"}, "QtyProteins": bson.M{"$size": "$Proteins"}, "_id": 0, "ScientificName": 1, "TaxId": 1}}}
	match := bson.D{{Key: "$match", Value: filter}}

	tt, err := db.table_basic.Aggregate(context.TODO(), mongo.Pipeline{match, groupStage})

	if err != nil {
		return orgTable, errorM.NewError("Database Error", errorM.StatusInternalServerError)
	}
	for tt.Next(context.TODO()) {
		var me str.OrganismTable
		var tax str.Taxonomy
		err = tt.Decode(&me)
		if err != nil {
			return orgTable, errorM.NewError("Can't Decode Result", errorM.StatusInternalServerError)
		}
		err = db.taxonomy_data.FindOne(context.TODO(), bson.D{{Key: "ScientificName", Value: me.ScientificName}}).Decode(&tax)
		if err != nil {
			return orgTable, errorM.NewError("Can't Find Taxon Id", errorM.StatusInternalServerError)
		}
		me.TaxId = tax.TaxID
		orgTable = append(orgTable, me)
	}

	return orgTable, nil
}

func (db *appDB) FindNucleotidesId(taxonId string) (str.NucleotideBasic, errorM.Errors) {
	var nBasic str.NucleotideBasic
	proj := options.FindOne().SetProjection(bson.D{{Key: "Proteins", Value: 0}})
	errM := db.table_basic.FindOne(context.TODO(), bson.D{{Key: "TaxId", Value: taxonId}}, proj).Decode(&nBasic)
	if errM != nil {
		return nBasic, errorM.NewError(errM.Error(), errorM.StatusBadRequest)
	}
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
