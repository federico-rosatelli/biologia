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
    <header>
        <hgroup>
            <h1>Micro Algae DB</h1>
        </hgroup>
        <nav>
            <ul class="right">
                <li>
                    <RouterLink to="/">
                        Home
					</RouterLink>
                </li>
                <li><a class="button-a" target=”_blank” href="https://github.com/federico-rosatelli/biologia">Github Project</a></li>
                <li><a class="button-a" target=”_blank” href="https://drive.google.com/drive/folders/19RXRHEb-7-O9gaUjXz5ho-Q2_HsbKlEW">Archive</a></li>
            </ul>
            <ul class="left">
                <li>
                    <RouterLink to="/taxonomy">
                        Taxonomy
					</RouterLink>
                </li>
                <li>
                    <RouterLink to="/taxon">
                        Search
					</RouterLink>
                </li>
            </ul>
        </nav>
    </header>
    <RouterView />
    <br/>
    <br/>
    <footer>
        <div id="banner" style="overflow:hidden; display: inline-block;">
            <div class="image">
                <img src="/public/images/Stemma_sapienza.png" style="width:30%; height:30%">
            </div>
        </div>
        <p>Authors: Team di Biologia Computazionale 2023, Università La Sapienza<br>
            <a href="mailto:">bio@bio.it</a>
        </p>
    </footer>
</template>

<style>
</style>
