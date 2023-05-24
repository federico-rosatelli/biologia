<!--IDEA: Vue delle ricerche con funzine che prende tipo di ricerca e input e ritorna tabella -->
<script>
import LoadingSpinner from '../components/LoadingSpinner.vue';


export default {
    data: function () {
        return {
            errormsg: null,
            loading: false,
            search: null,
            product: null,
            location: null,
            type: "scientific_name",
            response: []
        };
    },
    methods: {
        // ricerca con "stringa input" and "type": taxomID or scientific name
        async nucleotideSearch() {
            this.loading = true;
            this.errormsg = null;
            try {
                //  
                let response = await this.$axios.get("/taxon_term", {
                    params: {
                        search: this.search,
                        type: this.type,
                        product: this.product,
                        location: this.location,
                    }
                });
                this.response = response.data;
            }
            catch (e) {
                this.loading = false;
                this.errormsg = e.response.data;
                this.response = []
            }
            this.loading = false;
        },
        mouseOver(id,item){
          if (!item || item.length == 0){
            return
          }
          let fin = document.getElementById(id)
          console.log(fin);
          fin.style.visibility = "visible"
        },
        mouseNotOver(id,item){
          if (!item || item.length == 0){
            return
          }
          let fin = document.getElementById(id)
          console.log(fin);
          fin.style.visibility = "hidden"
        },
    },
    components: { LoadingSpinner }
}


</script>

<template>
  <div>
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
      <input type="text" id="search-input" placeholder="Input here" v-model="search" v-on:keyup.enter="nucleotideSearch()"/>
      <input type="text" id="search-input" placeholder="Filter for product" v-model="product" v-on:keyup.enter="nucleotideSearch()"/>
      <input type="text" id="search-input" placeholder="Filter for location" v-model="location" v-on:keyup.enter="nucleotideSearch()"/>
      <button type="button" @click="nucleotideSearch()"> Search</button>
    </div>
  </div>


<!--TABELLA -->
<br/>
<br/>
<br/>
<br/>
  <LoadingSpinner :loading="this.loading"/>
  <ErrorMsg v-if="errormsg" :msg="errormsg"></ErrorMsg>
  <div v-if="response.length > 0">
    <div>
      <h1>{{ response.length }} Results</h1>
    </div>
    <table id="table">
      <thead>
        <tr>
          <th>Scientific Name</th>
          <th>Qty Nucleotides</th>
          <th>Qty Proteins</th>
          <th>Qty Products</th>
          <th>Genomes</th>
          <th>Annotations</th>
          <th>Transcriptome</th>
          <th>Sra Wgs</th>
          <th>Sra Tran</th>
        </tr>
      </thead>
        <tr v-for="item in response" :key="item.TaxId">
          <td>
            <RouterLink :to="'/organism/'+item.TaxId + '/nucleotides'">
              {{ item.ScientificName }}
					  </RouterLink>
          </td>

          <td>
            <RouterLink :to="'/organism/'+item.TaxId+'/nucleotides'">
              {{ item.QtyNucleotides >= 9999 ? 9999+"+" : item.QtyNucleotides}}
					  </RouterLink>
          </td>

          <td>
            <RouterLink :to="'/organism/'+item.TaxId + '/proteins'">
              {{ item.QtyProteins >= 9999 ? 9999+"+" : item.QtyProteins}}
            </RouterLink>
          </td>

          <td>
            <RouterLink :to="'/organism/'+item.TaxId + '/nucleotides'">
              {{ item.QtyProducts >= 9999 ? 9999+"+" : item.QtyProducts}}
            </RouterLink>
          </td>
            
          <td @mouseleave="mouseNotOver('genome-'+item.TaxId,item.Genomes)" @mouseover="mouseOver('genome-'+item.TaxId,item.Genomes)">
            {{ !item.Genomes || item.Genomes.length > 0 ? 'Link' : '-' }}
            <div class="dropdown-genome">
              <div :id="'genome-'+item.TaxId" class="dropdown-content-genome">
                <a v-if="item.Genomes && item.Genomes.length > 0" :href="item.Genomes[0].Link" target=”_blank”>
                  <div :style="item.Genomes[0].GBFF ? 'background-color:green;' : 'background-color:red;'">
                    GCFF
                  </div>
                  <div :style="item.Genomes[0].FNA ? 'background-color:green;' : 'background-color:red;'">
                    FNA
                  </div>
                  <div :style="item.Genomes[0].GFF ? 'background-color:green;' : 'background-color:red;'">
                    GFF
                  </div>
                </a>
              </div>
            </div>
            
          </td>
          <td>{{ item.Annotations }}</td>  <!--link download gff -->
          <td>{{ item.Trascriptome }}</td> <!--link tsa ncbi  -->
          <td>{{ item.SraWgs }}</td>
          <td>{{ item.SraTran }}</td>
        </tr>
    </table>
  </div>


</template>


<style>
.input {
  margin-bottom: 10px;
}

.dropdown {
  margin-bottom: 20px;
}

button {
  padding: 10px;
  background-color: #003399;
  color: white;
  border-radius: 10px;
  border: none;
  cursor: pointer;
}

button:hover {
  background-color: #0b11bb;
}


.row{
    display: flex;
    gap: 5px;
   }
   #table {
    border-collapse: collapse;
    width: 100%;
   }
   #table td, #table th {
    border: 1px solid #6cace1;
    padding: 20px;
    align-content: center;
   }
   #table .title{
     font-size: 25px;
     background-color: #999999;
     text-align: center;
   }
   /* #table tr:nth-child(even){background-color: #88a3f7;} */
   table tr:hover {background-color: rgb(134, 220, 239,.18);}
   #table th {
    padding-top: 12px;
    padding-bottom: 12px;
    text-align: left;
    background-color: #04aa6d;
    
   }

.input input[type="text"] {
    width: 40%;
    padding: 10px;
    left: 50%;
    border: 1px solid #dddddd;
    margin-bottom: 15px;
    box-sizing:border-box;
}

.dropdown-genome {
  position: relative;
  display: block;
  margin-right: 25px;
}

.dropdown-content-genome {
  visibility: hidden;
  position: absolute;
  min-width: 100px;
  box-shadow: 0px 8px 16px 0px rgba(0,0,0,0.2);
  z-index: 1;
  align-content: center;
}

.dropdown-genome:hover .dropdown-content-genome {
  display: block;
}
</style>