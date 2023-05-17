import {createRouter, createWebHistory} from 'vue-router'
import SearchView from '../views/SearchView.vue'

const router = createRouter({
	history: createWebHistory(),
	routes: [
		{path: '/', component: SearchView},
		{path: '/taxon', component: SearchView},

	]
})

export default router
