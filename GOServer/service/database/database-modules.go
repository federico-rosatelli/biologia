package database

import (
	"context"

	"go.mongodb.org/mongo-driver/bson"
)

type Prova1 struct {
	TaxId string
}

func (db *appDB) FindTaxon(search string, typeS string) (p Prova1, err error) {

	filter := bson.D{{Key: "ScientificName", Value: search}}
	err = db.taxonomy_data.FindOne(context.TODO(), filter).Decode(&p)

	return p, err
}
