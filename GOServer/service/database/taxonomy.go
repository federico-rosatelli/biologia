package database

import (
	str "GOServer/service/structures"
	"context"

	"go.mongodb.org/mongo-driver/bson"
)

func (db *appDB) FindTaxon(search string, typeS string) (str.Taxonomy, error) {
	var t str.Taxonomy
	filter := bson.M{"ScientificName": search}
	err := db.taxonomy_data.FindOne(context.TODO(), filter).Decode(&t)
	return t, err
}

func (db *appDB) FindTaxonTree(search string) (str.TaxonomyTree, error) {
	var t str.TaxonomyTree
	filter := bson.M{"ScientificName": search}
	err := db.taxonomy_tree.FindOne(context.TODO(), filter).Decode(&t)
	return t, err
}
