package api

import (
	"GOServer/service/ErrManager"
	str "GOServer/service/structures"
)

func (rt *_router) GetAnalysis() ([]str.SraBioProject, ErrManager.Errors) {
	var an []str.SraBioProject
	nProj, err := rt.db.FindAnalysis()
	if err != nil {
		return an, err
	}
	var sraProjects []str.SraBioProject
	srr := make(map[string]str.SraBioProject)
	for _, pro := range nProj {
		projId := pro["BioProject"].(string)
		val, ok := srr[projId]
		if !ok {

			sra := str.SraBioProject{
				TaxId: pro["TaxId"].(string),

				BioProjectId: pro["BioProject"].(string),
				BioSamples: []struct {
					BioSampleId string
					Experiments []struct {
						ExperimentId string
						Runs         []string
					}
				}{},
			}
			if _, ok := pro["ScientificName"].(string); ok {
				sra.ScientificName = pro["ScientificName"].(string)
			}
			if _, ok := pro["Abstract"].(string); ok {
				sra.Comment = pro["Abstract"].(string)
			}
			sra.BioSamples = append(sra.BioSamples, struct {
				BioSampleId string
				Experiments []struct {
					ExperimentId string
					Runs         []string
				}
			}{
				BioSampleId: pro["BioSample"].(string),
			})
			sra.BioSamples[0].Experiments = append(sra.BioSamples[0].Experiments, struct {
				ExperimentId string
				Runs         []string
			}{ExperimentId: pro["ExperimentCode"].(string)})
			sra.BioSamples[0].Experiments[0].Runs = append(sra.BioSamples[0].Experiments[0].Runs, pro["Run"].(string))
			srr[projId] = sra
		} else {
			isInsideSample := false

			for k, sample := range val.BioSamples {
				if sample.BioSampleId == pro["BioSample"].(string) {
					isInsideSample = true
					isInsideExample := false
					for j, expreiment := range sample.Experiments {
						if expreiment.ExperimentId == pro["ExperimentCode"].(string) {
							isInsideExample = true
							srr[projId].BioSamples[k].Experiments[j].Runs = append(srr[projId].BioSamples[k].Experiments[j].Runs, pro["Run"].(string))
						}
					}
					if !isInsideExample {
						srr[projId].BioSamples[k].Experiments = append(srr[projId].BioSamples[k].Experiments, struct {
							ExperimentId string
							Runs         []string
						}{
							ExperimentId: pro["ExperimentCode"].(string),
						})
						srr[projId].BioSamples[k].Experiments[len(srr[projId].BioSamples[k].Experiments)-1].Runs = append(srr[projId].BioSamples[k].Experiments[len(srr[projId].BioSamples[k].Experiments)-1].Runs, pro["Run"].(string))
					}

				}
			}
			if !isInsideSample {
				prov := struct {
					BioSampleId string
					Experiments []struct {
						ExperimentId string
						Runs         []string
					}
				}{BioSampleId: pro["BioSample"].(string)}

				prov.Experiments = append(prov.Experiments, struct {
					ExperimentId string
					Runs         []string
				}{ExperimentId: pro["ExperimentCode"].(string)})
				prov.Experiments[0].Runs = append(prov.Experiments[0].Runs, pro["Run"].(string))
				value := srr[projId]
				value.BioSamples = append(value.BioSamples, prov)
				srr[projId] = value
			}

		}
	}
	for _, key := range srr {
		sraProjects = append(sraProjects, key)
	}
	return sraProjects, nil
}

// func (a myAnalysis) mutateStruct() (str.SraBioProject, ErrManager.Errors) {

// }
