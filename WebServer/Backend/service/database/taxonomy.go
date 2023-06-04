package database

import (
	errorM "GOServer/service/ErrManager"
	str "GOServer/service/structures"
	"context"

	"go.mongodb.org/mongo-driver/bson"
)

func (db *appDB) FindTaxon(search string) (str.Taxonomy, errorM.Errors) {
	var t str.Taxonomy
	filter := bson.M{"ScientificName": search}
	err := db.taxonomy_data.FindOne(context.TODO(), filter).Decode(&t)
	if err != nil {
		return t, errorM.NewError("Can't find Taxon info", errorM.StatusBadRequest)
	}
	return t, nil
}

func (db *appDB) FindTaxonTree(search string, typeN string) (str.TaxonomyTree, errorM.Errors) {
	var t str.TaxonomyTree
	filter := bson.M{}
	if typeN == "id" {
		filter["TaxId"] = search
	} else if typeN == "scientific_name" {
		filter["ScientificName"] = bson.D{{Key: "$regex", Value: search}, {Key: "$options", Value: "i"}}
	} else {
		return t, errorM.NewError("Only id|scientific_name allowed", errorM.StatusBadRequest)
	}
	err := db.taxonomy_tree.FindOne(context.TODO(), filter).Decode(&t)
	if err != nil {
		return t, errorM.NewError("Can't find Taxon info", errorM.StatusBadRequest)
	}
	for i := 0; i < len(t.SubClasses); i++ {
		var t1 str.TaxonomyTree
		filter := bson.M{"TaxId": t.SubClasses[i].TaxId}
		_ = db.taxonomy_tree.FindOne(context.TODO(), filter).Decode(&t1)
		for k := 0; k < len(t1.SubClasses); k++ {
			var t2 str.TaxonomyTree
			filter := bson.M{"TaxId": t1.SubClasses[k].TaxId}
			_ = db.taxonomy_tree.FindOne(context.TODO(), filter).Decode(&t2)
			for j := 0; j < len(t2.SubClasses); j++ {
				var t3 str.TaxonomyTree
				filter := bson.M{"TaxId": t2.SubClasses[j].TaxId}
				_ = db.taxonomy_tree.FindOne(context.TODO(), filter).Decode(&t3)

				t2.SubClasses[j] = t3
			}
			t1.SubClasses[k] = t2
		}
		t.SubClasses[i] = t1
	}
	return t, nil
}
