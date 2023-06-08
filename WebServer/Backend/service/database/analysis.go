package database

import (
	errorM "GOServer/service/ErrManager"
	"context"

	"go.mongodb.org/mongo-driver/bson"
	"go.mongodb.org/mongo-driver/mongo"
)

// filter = {
// 	"STUDY.center_name": "BioProject",
// 	"SAMPLE.IDENTIFIERS.EXTERNAL_ID.namespace": "BioSample",
// }

// proj = {
// 	"_id":0,
// 	"TaxId":"$SAMPLE.SAMPLE_NAME.TAXON_ID.value",
// 	"BioProject":"$STUDY.alias",
// 	"BioSample":"$SAMPLE.IDENTIFIERS.EXTERNAL_ID.value",
// 	"ExperimentCode":"$EXPERIMENT.accession",
// 	"SRR":"$RUN_SET.RUN.accession",
// 	"Abstract":"$STUDY.DESCRIPTOR.STUDY_ABSTRACT.value",
// 	"LinkSRA":"$RUN_SET.RUN.SRAFiles"
// }

func (db *appDB) FindAnalysis() ([]map[string]interface{}, errorM.Errors) {
	filter := bson.D{{
		Key: "$match", Value: bson.M{
			"STUDY.center_name":                        "BioProject",
			"SAMPLE.IDENTIFIERS.EXTERNAL_ID.namespace": "BioSample",
		},
	}}
	proj := bson.D{{
		Key: "$project", Value: bson.M{
			"_id":            0,
			"ScientificName": "$SAMPLE.SAMPLE_NAME.SCIENTIFIC_NAME.value",
			"TaxId":          "$SAMPLE.SAMPLE_NAME.TAXON_ID.value",
			"BioProject":     "$STUDY.alias",
			"BioSample":      "$SAMPLE.IDENTIFIERS.EXTERNAL_ID.value",
			"ExperimentCode": "$EXPERIMENT.accession",
			"Run":            "$RUN_SET.RUN.accession",
			"Abstract":       "$STUDY.DESCRIPTOR.STUDY_ABSTRACT.value",
		},
	}}
	dataSeq, errM := db.sequences_data.Aggregate(context.TODO(), mongo.Pipeline{filter, proj})
	if errM != nil {
		return nil, errorM.NewError(errM.Error(), errorM.StatusInternalServerError)
	}
	var seq_data []map[string]interface{}
	for dataSeq.Next(context.TODO()) {
		var seq map[string]interface{}
		err := dataSeq.Decode(&seq)
		if err != nil {
			return nil, errorM.NewError("Can't Decode Result", errorM.StatusInternalServerError)
		}
		seq_data = append(seq_data, seq)
	}
	return seq_data, nil
}
