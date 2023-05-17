import {createApp, reactive} from 'vue'
// import { BModal, BButton } from 'bootstrap-vue-3'

import App from './App.vue'
import router from './router'
import axios from './services/axios.js';
import ErrorMsg from './components/ErrorMsg.vue'
import './assets/dashboard.css'
import './assets/main.css'
import './assets/login.css'
import './assets/modal.css'
import './assets/dropdown.css'
console.log("PROVAAA");
const app = createApp(App)
app.config.globalProperties.$axios = axios;
app.component("ErrorMsg", ErrorMsg);
app.use(router)
//app.use(BModal)
// app.component("b-modal",BModal)
// app.component("b-button",BButton)
//app.use(IconsPlugin)

//app.component('b-modal', BModal)
app.mount('#app')

