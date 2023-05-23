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

	FindTaxon(search string) (str.Taxonomy, error)
	FindTaxonTree(taxId string) (str.TaxonomyTree, error)
	TableOrganism(search string, typeS string, listTaxId []string) ([]str.OrganismTable, errorM.Errors)
	TableOrganismByProduct(product string) ([]string, errorM.Errors)
	TableOrganismByLocation(location string) ([]string, errorM.Errors)
	FindNucleotidesId(taxonId string) (str.TableBasic, errorM.Errors)
	FindNucleotideTableByLocus(taxId string, locus string) (str.TableComplete, errorM.Errors)
	FindNucleotideByLocus(locus string) (str.Nucleotide, errorM.Errors)
	FindProteinsId(taxonId string) (str.TableBasic, errorM.Errors)
	FindProteinTableByLocus(taxId string, locus string) (str.TableComplete, errorM.Errors)
	FindProteinByLocus(locus string) (str.Protein, errorM.Errors)
}

type appDB struct {
	client          *mongo.Client
	nucleotide_data *mongo.Collection
	protein_data    *mongo.Collection
	table_basic     *mongo.Collection
	table_complete  *mongo.Collection
	taxonomy_data   *mongo.Collection
	taxonomy_tree   *mongo.Collection
}

func InitDatabase(client *mongo.Client) (AppDatabase, error) {
	collectionNucleotide := client.Database("Biologia").Collection("nucleotide_data")
	if collectionNucleotide == nil {
		return nil, errors.New("error Creating users Collection")
	}
	collectionProtein := client.Database("Biologia").Collection("protein_data")
	if collectionProtein == nil {
		return nil, errors.New("error Creating users Collection")
	}
	collectionTableBasic := client.Database("Biologia").Collection("table_basic")
	if collectionTableBasic == nil {
		return nil, errors.New("error Creating users Collection")
	}
	collectionTableComplete := client.Database("Biologia").Collection("table_complete")
	if collectionTableComplete == nil {
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
		protein_data:    collectionProtein,
		table_basic:     collectionTableBasic,
		table_complete:  collectionTableComplete,
		taxonomy_data:   collectionTaxonomy,
		taxonomy_tree:   collectionTaxonomyTree,
	}, nil
}
