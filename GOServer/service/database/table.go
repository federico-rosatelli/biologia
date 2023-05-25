package database

import (
	errorM "GOServer/service/ErrManager"
	str "GOServer/service/structures"
	"context"
	"strings"

	"go.mongodb.org/mongo-driver/bson"
	"go.mongodb.org/mongo-driver/mongo"
	"go.mongodb.org/mongo-driver/mongo/options"
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

func (db *appDB) WelcomeMarkDown() (str.Markdown, errorM.Errors) {
	var mark str.Markdown
	err := db.markdown.FindOne(context.TODO(), bson.D{{}}).Decode(&mark)
	if err != nil {
		return mark, errorM.NewError("Can't find Markdown info", errorM.StatusBadRequest)
	}
	return mark, nil
}

func (db *appDB) TableOrganismByProduct(product string) ([]string, errorM.Errors) {
	var listTaxId []string
	proj := options.Find().SetProjection(bson.D{{Key: "TaxId", Value: 1}})
	finder := bson.M{"Products.ProductName": bson.M{
		"$regex": product, "$options": "i"}}
	cursor, errM := db.table_basic.Find(context.TODO(), finder, proj)
	if errM != nil {
		return listTaxId, errorM.NewError(errM.Error(), errorM.StatusBadRequest)
	}

	for cursor.Next(context.TODO()) {
		var tabl str.TableBasic
		err := cursor.Decode(&tabl)
		if err != nil {
			return listTaxId, errorM.NewError("Can't Decode Result", errorM.StatusInternalServerError)
		}
		listTaxId = append(listTaxId, tabl.TaxId)
	}

	return listTaxId, nil
}

func (db *appDB) TableOrganismByLocation(location string) ([]string, errorM.Errors) {
	var listTaxId []string
	proj := options.Find().SetProjection(bson.D{{Key: "GBSeq_feature-table", Value: 1}})
	finder := bson.M{"GBSeq_feature-table.GBFeature_key": "source",
		"GBSeq_feature-table.GBFeature_quals.GBQualifier_name": "country",
		"GBSeq_feature-table.GBFeature_quals.GBQualifier_value": bson.M{
			"$regex": location, "$options": "i"}}
	cursor, errM := db.nucleotide_data.Find(context.TODO(), finder, proj)
	if errM != nil {
		return listTaxId, errorM.NewError(errM.Error(), errorM.StatusBadRequest)
	}

	for cursor.Next(context.TODO()) {
		var nucleo str.Nucleotide
		err := cursor.Decode(&nucleo)
		if err != nil {
			return listTaxId, errorM.NewError("Can't Decode Result", errorM.StatusInternalServerError)
		}
		taxId := ""
		quals := nucleo.GBSeqFeatureTable[0].GBFeatureQuals
		for i := 0; i < len(quals); i++ {
			if quals[i].GBQualifierName == "db_xref" {
				taxId = strings.Split(quals[i].GBQualifierValue, ":")[1]
			}
		}
		if taxId != "" {
			listTaxId = append(listTaxId, taxId)
		}
	}

	return listTaxId, nil
}
