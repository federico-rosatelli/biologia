import {createRouter, createWebHistory} from 'vue-router'
import SearchView from '../views/search.vue'
import HomeView from '../views/Home.vue'

const router = createRouter({
	history: createWebHistory(),
	routes: [
		{path: '/', component: HomeView},
		{path: '/taxon', component: SearchView},

	]
})

export default router
