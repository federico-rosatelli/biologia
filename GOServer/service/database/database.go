package database

import (
	str "GOServer/service/structures"
	"errors"

	"go.mongodb.org/mongo-driver/mongo"
)

type AppDatabase interface {
	// GetName() (string, error)
	// SetName(name string) error

	FindTaxon(search string, typeS string) (str.Taxonomy, error)
	FindTaxonTree(search string) (str.TaxonomyTree, error)
	FindNucleotide(search string, typeS string) ([]map[string]interface{}, error)
}

type appDB struct {
	client          *mongo.Client
	nucleotide_data *mongo.Collection
	taxonomy_data   *mongo.Collection
	taxonomy_tree   *mongo.Collection
}

func InitDatabase(client *mongo.Client) (AppDatabase, error) {
	collectionNucleotide := client.Database("Biologia").Collection("nucleotide_data")
	if collectionNucleotide == nil {
		return nil, errors.New("error Creating users Collection")
	}
	collectionTaxonomy := client.Database("Biologia").Collection("taxonomy_data")
	if collectionTaxonomy == nil {
		return nil, errors.New("error Creating users Collection")
	}
	collectionTaxonomyTree := client.Database("Biologia").Collection("taxonomy_tree")
	if collectionTaxonomyTree == nil {
		return nil, errors.New("error Creating users Collection")
	}
	return &appDB{
		client:          client,
		nucleotide_data: collectionNucleotide,
		taxonomy_data:   collectionTaxonomy,
		taxonomy_tree:   collectionTaxonomyTree,
	}, nil
}
