import {createRouter, createWebHistory} from 'vue-router'
import SearchView from '../views/search.vue'
import HomeView from '../views/Home.vue'
import OrganismView from '../views/OrganismView.vue'

const router = createRouter({
	history: createWebHistory(),
	routes: [
		{path: '/', component: HomeView},
		{path: '/taxon', component: SearchView},
		{path: '/organism/:taxid/nucleotides', component: OrganismView},
		{path: '/organism/:taxid/proteins', component: OrganismView},

	]
})

export default router
