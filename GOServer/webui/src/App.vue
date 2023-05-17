<script>
import { RouterLink, RouterView } from 'vue-router'
export default{
    data:()=>{
        return {
            token: null,
            some_data: null,
            errormsg: null,
			dataBool: false,
        }
    },
    methods:{
        async info(){
			// let response = await this.$axios.get("/");
			// this.some_data = response.data;

        },
    },
    async mounted(){
        await this.info()
    }

}
</script>

<template>

	
	<div class="container-fluid">
		<div class="row">
			<main>
				<RouterView />
			</main>
		</div>
	</div>
</template>

<style>
</style>
