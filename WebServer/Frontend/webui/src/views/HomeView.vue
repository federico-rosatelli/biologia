<script>
import { marked } from 'marked';
export default {
	data: function() {
		return {
			errormsg: null,
			loading: false,
			some_data: ""
		}
	},
	methods: {
        data(){
        },
	},
	computed:{
		compiledMarkdown() {
			
			return marked(this.some_data);
		}
	},
	async mounted(){
		try{
            let response = await this.$axios.get(`/`);
            this.some_data=response.data.Text
        } catch(e){
            this.errormsg = e.response.data
        }
	}
}
</script>

<template>
	<div v-html="compiledMarkdown"></div>
</template>