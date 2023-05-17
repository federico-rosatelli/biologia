package database

import (
	"errors"

	"go.mongodb.org/mongo-driver/mongo"
)

type AppDatabase interface {
	// GetName() (string, error)
	// SetName(name string) error

	FindTaxon(search string, typeS string) (p Prova1, err error)
}

type appDB struct {
	client          *mongo.Client
	nucleotide_data *mongo.Collection
	taxonomy_data   *mongo.Collection
}

func InitDatabase(client *mongo.Client) (AppDatabase, error) {
	collectionNucleotide := client.Database("biologia").Collection("nucleotide_data")
	if collectionNucleotide == nil {
		return nil, errors.New("error Creating users Collection")
	}
	collectionTaxonomy := client.Database("biologia").Collection("taxonomy_data")
	if collectionTaxonomy == nil {
		return nil, errors.New("error Creating users Collection")
	}
	return &appDB{
		client:          client,
		nucleotide_data: collectionNucleotide,
		taxonomy_data:   collectionTaxonomy,
	}, nil
}
