<!--IDEA: Vue delle ricerche con funzine che prende tipo di ricerca e input e ritorna tabella -->
<script>

export default {
	data: function() {
		return {
			errormsg: null,
			loading: false,
            search: null,
            data: null,
            type: "scientific_name",

		}
	},
    
	methods: {
        async searchTax(){
            this.errormsg = null
			this.loading = true
            try{
                let data = await this.$axios.get(`taxonomy_tree`, {
                    params: {
                        search:this.search,
                        type: this.type
                    }
                });
                this.data = data.data
            } catch(e){
                this.errormsg = e.response.data;
                this.loading = false
            }
            this.loading = false
        },
    },
    async mounted(){
        let taxId = this.$route.params.taxid;
        console.log(taxId);
        if (!taxId){
            return
        }
        this.loading = true
        try{
            let data = await this.$axios.get(`taxonomy_tree`, {
                params: {
                        search:taxId,
                        type: 'id'
                }
            });
            this.data = data.data
        } catch(e){
            this.errormsg = e.response.data;
            this.loading = false
        }
        this.loading = false
    },
}


</script>

<template>
    <LoadingSpinner :loading="this.loading"/>
    <div class="dropdown">
      <label for="dropdown">
        <h1>
          Type:
        </h1>
      </label>
      <select id="dropdown" name="drop minimal" v-model="type">
        <option value="scientific_name" selected>Scientific Name</option>
        <option value="id">Taxon Id</option>
      </select>
    </div>
    <div class="input">
      <input type="text" id="search-input" placeholder="Input here" v-model="search" v-on:keyup.enter="searchTax()"/>
      <button type="button" @click="searchTax()"> Search</button>
    </div>
    <ErrorMsg v-if="errormsg" :msg="errormsg"></ErrorMsg>
    <div v-if="data">
        <h2 class="scientificname">
            {{data.ScientificName}}
        </h2>
        <br/>
        <ul v-for="t1 in data.SubClasses">
            <a v-if="t1.Rank == 'species' || t1.Rank == 'varietas'" :href="'/organism/'+t1.TaxId+'/nucleotides'">
                {{ t1.ScientificName }} {{ t1.Rank }}
            </a>
            <a v-else :href="'/taxonomy/'+t1.TaxId">
                {{ t1.ScientificName }} {{ t1.Rank }}
            </a>
            <ul v-for="t2 in t1.SubClasses">
                <a v-if="t2.Rank == 'species' || t2.Rank == 'varietas'" :href="'/organism/'+t2.TaxId+'/nucleotides'">
                {{ t2.ScientificName }} {{ t2.Rank }}
                </a>
                <a v-else :href="'/taxonomy/'+t2.TaxId">
                    {{ t2.ScientificName }} {{ t2.Rank }}
                </a>
                <ul v-for="t3 in t2.SubClasses">
                    <a v-if="t3.Rank == 'species' || t3.Rank == 'varietas'" :href="'/organism/'+t3.TaxId+'/nucleotides'">
                        {{ t3.ScientificName }} {{ t3.Rank }}
                    </a>
                    <a v-else :href="'/taxonomy/'+t3.TaxId">
                        {{ t3.ScientificName }} {{ t3.Rank }}
                    </a>
                    <ul v-for="t4 in t3.SubClasses">
                        <a v-if="t4.Rank == 'species' || t4.Rank == 'varietas'" :href="'/organism/'+t4.TaxId+'/nucleotides'">
                            {{ t4.ScientificName }} {{ t4.Rank }}
                        </a>
                        <a v-else :href="'/taxonomy/'+t4.TaxId">
                            {{ t4.ScientificName }} {{ t4.Rank }}
                        </a>
                    </ul>
                </ul>
            </ul>
        </ul>
    </div>

</template>

<style>
.scientificname{
    position: absolute;
    text-align: left;
}
</style>
