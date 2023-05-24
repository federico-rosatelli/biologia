package database

import (
	errorM "GOServer/service/ErrManager"
	str "GOServer/service/structures"
	"context"

	"go.mongodb.org/mongo-driver/bson"
	"go.mongodb.org/mongo-driver/mongo"
)

func (db *appDB) TableOrganism(search string, typeS string, listTaxId []string) ([]str.OrganismTable, errorM.Errors) {
	var orgTable []str.OrganismTable
	filter := bson.M{}
	if typeS == "id" {
		filter["TaxId"] = search
	} else if typeS == "scientific_name" {
		filter["ScientificName"] = bson.D{{Key: "$regex", Value: search}, {Key: "$options", Value: "i"}}
	} else {
		return orgTable, errorM.NewError("Only id|scientific_name allowed", errorM.StatusBadRequest)
	}
	if len(listTaxId) > 0 {
		filter["TaxId"] = bson.D{{Key: "$in", Value: listTaxId}}
	}

	groupStage := bson.D{{Key: "$project", Value: bson.M{"QtyNucleotides": bson.M{"$size": "$Nucleotides"}, "QtyProteins": bson.M{"$size": "$Proteins"}, "QtyProducts": bson.M{"$size": "$Products"}, "_id": 0, "ScientificName": 1, "TaxId": 1, "Genomes": 1}}}
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
