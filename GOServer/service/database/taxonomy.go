package database

import (
	str "GOServer/service/structures"
	"context"

	"go.mongodb.org/mongo-driver/bson"
)

func (db *appDB) FindTaxon(search string) (str.Taxonomy, error) {
	var t str.Taxonomy
	filter := bson.M{"ScientificName": search}
	err := db.taxonomy_data.FindOne(context.TODO(), filter).Decode(&t)
	return t, err
}

func (db *appDB) FindTaxonTree(taxId string) (str.TaxonomyTree, error) {
	var t str.TaxonomyTree
	filter := bson.M{"TaxId": taxId}
	err := db.taxonomy_tree.FindOne(context.TODO(), filter).Decode(&t)
	if err != nil {
		return t, err
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
	return t, err
}
