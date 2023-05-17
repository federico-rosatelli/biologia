package database

import (
	"context"
	"errors"

	"go.mongodb.org/mongo-driver/bson"
)

func (db *appDB) FindNucleotide(search string, typeS string) ([]map[string]interface{}, error) {
	var ttt []map[string]interface{}
	filter := bson.M{}
	if typeS == "id" {
		filter["GBSeq_locus"] = search
	} else if typeS == "scientific_name" {
		filter["GBSeq_source"] = search
	} else {
		return ttt, errors.New("only id|scientific_name allowed")
	}
	// db.collection.aggregate( [
	// 	{ $group: { _id: null, myCount: { $sum: 1 } } },
	// 	{ $project: { _id: 0 } }
	//  ] )
	// { "$group": {
	// 	"_id": "$_id",
	// 	"name": { "$first": "$name" },  //$first accumulator
	// 	"count": { "$sum": 1 },  //$sum accumulator
	// 	"totalValue": { "$sum": "$value" }  //$sum accumulator
	//   }}
	// groupStage := bson.D{{Key: "$group", Value: bson.M{"_id": "$GBSeq_locus", "myCount": bson.M{"$sum": 1}}}}
	// finder := bson.D{{Key: "$match", Value: bson.D{{Key: "GBSeq_locus", Value: search}}}}
	// //err := db.nucleotide_data.FindOne(context.TODO(), filter).Decode(&nucleo)
	tt, err := db.nucleotide_data.Find(context.TODO(), bson.D{{Key: "GBSeq_locus", Value: search}})
	// tt, err := db.nucleotide_data.Aggregate(context.TODO(), mongo.Pipeline{groupStage, finder})
	print("MA QUII??")
	if err != nil {
		print(err.Error(), "LO SAPEVO")
		return ttt, err
	}
	print("EH MA QUA ARRIVA")
	if err != nil {
		print(err.Error(), "  LO SAPEVO123")
		return ttt, err
	}
	for tt.Next(context.TODO()) {
		var me map[string]interface{}
		err = tt.Decode(&me)
		if err != nil {
			print(err.Error())
			return ttt, err
		}
		ttt = append(ttt, me)
	}
	print(ttt)
	return ttt, err
}
