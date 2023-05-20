import {createRouter, createWebHistory} from 'vue-router'
import SearchView from '../views/search.vue'
import HomeView from '../views/Home.vue'
import NucleosView from '../views/NucleosView.vue'

const router = createRouter({
	history: createWebHistory(),
	routes: [
		{path: '/', component: HomeView},
		{path: '/taxon', component: SearchView},
		{path: '/organism/:taxid/nucleotides', component: NucleosView},

	]
})

export default router
