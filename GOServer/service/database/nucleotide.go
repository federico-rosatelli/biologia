package database

import (
	errorM "GOServer/service/ErrManager"
	str "GOServer/service/structures"
	"context"

	"go.mongodb.org/mongo-driver/bson"
	"go.mongodb.org/mongo-driver/mongo"
)

func (db *appDB) TableOrganism(search string, typeS string) ([]str.OrganismTable, errorM.Errors) {
	var orgTable []str.OrganismTable
	filter := bson.M{}
	if typeS == "id" {
		filter["GBSeq_locus"] = bson.D{{Key: "$elemMatch", Value: bson.D{{Key: "GBSeq_locus", Value: search}}}}
	} else if typeS == "scientific_name" {
		filter["ScientificName"] = bson.D{{Key: "$regex", Value: search}}
	} else {
		return orgTable, errorM.NewError("Only id|scientific_name allowed", errorM.StatusBadRequest)
	}

	groupStage := bson.D{{Key: "$project", Value: bson.M{"QtyNucleotides": bson.M{"$size": "$GBSeq_locus"}, "_id": 0, "ScientificName": 1}}}
	match := bson.D{{Key: "$match", Value: filter}}

	tt, err := db.nucleotide_base.Aggregate(context.TODO(), mongo.Pipeline{match, groupStage})

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
	var tax str.Taxonomy
	errM := db.taxonomy_data.FindOne(context.TODO(), bson.D{{Key: "TaxId", Value: taxonId}}).Decode(&tax)
	if errM != nil {
		return nBasic, errorM.NewError(errM.Error(), errorM.StatusBadRequest)
	}
	errM = db.nucleotide_base.FindOne(context.TODO(), bson.D{{Key: "ScientificName", Value: tax.ScientificName}}).Decode(&nBasic)
	if errM != nil {
		return nBasic, errorM.NewError(errM.Error(), errorM.StatusBadRequest)
	}
	return nBasic, nil

}
