package database

import (
	"context"

	"go.mongodb.org/mongo-driver/mongo/readpref"
)

func (db *appDB) Ping() (err error) {
	err = db.client.Ping(context.TODO(), readpref.Primary())
	return err
}
