package database

import (
	errorM "GOServer/service/ErrManager"
	str "GOServer/service/structures"
	"errors"

	"go.mongodb.org/mongo-driver/mongo"
)

type AppDatabase interface {
	// GetName() (string, error)
	// SetName(name string) error

	FindTaxon(search string, typeS string) (str.Taxonomy, error)
	FindTaxonTree(search string) (str.TaxonomyTree, error)
	TableOrganism(search string, typeS string) ([]str.OrganismTable, errorM.Errors)
	FindNucleotidesId(taxonId string) (str.NucleotideBasic, errorM.Errors)
}

type appDB struct {
	client          *mongo.Client
	nucleotide_data *mongo.Collection
	table_basic     *mongo.Collection
	taxonomy_data   *mongo.Collection
	taxonomy_tree   *mongo.Collection
}

func InitDatabase(client *mongo.Client) (AppDatabase, error) {
	collectionNucleotide := client.Database("Biologia").Collection("nucleotide_data")
	if collectionNucleotide == nil {
		return nil, errors.New("error Creating users Collection")
	}
	collectionTableBasic := client.Database("Biologia").Collection("table_basic")
	if collectionTableBasic == nil {
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
		table_basic:     collectionTableBasic,
		taxonomy_data:   collectionTaxonomy,
		taxonomy_tree:   collectionTaxonomyTree,
	}, nil
}
