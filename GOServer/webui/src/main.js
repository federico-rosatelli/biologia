import {createApp, reactive} from 'vue'
// import { BModal, BButton } from 'bootstrap-vue-3'

import App from './App.vue'
import router from './router'
import axios from './services/axios.js';
import ErrorMsg from './components/ErrorMsg.vue'
import NucleoComp from './components/NucleotideComp.vue'
import LoadingSpinner from './components/LoadingSpinner.vue'
import './assets/dashboard.css'
import './assets/main.css'
import './assets/login.css'
import './assets/modal.css'
import './assets/dropdown.css'
import './assets/monitor.css'
console.log("PROVAAA");
const app = createApp(App)
app.config.globalProperties.$axios = axios;
app.component("ErrorMsg", ErrorMsg);
app.component("NucleoComp", NucleoComp);
app.component("LoadingSpinner", LoadingSpinner);
app.use(router)
app.mount('#app')

