package api

import (
	"GOServer/service/api/ErrManager"
	str "GOServer/service/structures"
)

func (rt *_router) GetOrganism(query string, nType string) (str.OrganismTable, ErrManager.Errors) {
	var o str.OrganismTable
	return o, nil
}
