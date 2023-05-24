import {createRouter, createWebHistory} from 'vue-router'
import SearchTable from '../views/SearchTable.vue'
import HomeView from '../views/Home.vue'
import Taxonomy from '../views/TaxonomyView.vue'
import OrganismView from '../views/OrganismView.vue'

const router = createRouter({
	history: createWebHistory(),
	routes: [
		{path: '/', component: HomeView},
		{path: '/taxonomy', component: Taxonomy},
		{path: '/taxonomy/:taxid', component: Taxonomy},
		{path: '/taxon', component: SearchTable},
		{path: '/organism/:taxid/nucleotides', component: OrganismView},
		{path: '/organism/:taxid/proteins', component: OrganismView},

	]
})

export default router
