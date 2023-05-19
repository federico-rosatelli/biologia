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
            <h1>Danio Rerio DB</h1>
            <h2>Powered by MongoDB</h2>
        </hgroup>
        <nav>
            <ul class="right">
                <li><a href="/">Home</a></li>
                <li><a class="button-a" target=”_blank” href="https://github.com/federico-rosatelli/biologia">Github Project</a></li>
                <li><a class="button-a" target=”_blank” href="https://drive.google.com/drive/folders/19RXRHEb-7-O9gaUjXz5ho-Q2_HsbKlEW">Archive</a></li>
            </ul>
            <ul class="left">
                <li><a href="/taxon">Taxon</a></li>
            </ul>
        </nav>
    </header>
    <RouterView />
    <footer>
    <div id="banner" style="overflow:hidden; display: inline-block;">
      <div class="image">
        <img src="images/Stemma_sapienza.png" style="width:30%; height:30%">
      </div>
    </div>
    <p>Authors: Team di Biologia Computazionale 2023, Università La Sapienza<br>
    <a href="mailto:">bio@bio.it</a></p>
  </footer>
</template>

<style>
</style>
