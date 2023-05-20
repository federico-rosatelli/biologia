<!--IDEA: Vue delle ricerche con funzine che prende tipo di ricerca e input e ritorna tabella -->
<script>

export default {
	data: function() {
		return {
			errormsg: null,
			loading: false,
            search: null,
            type: null,
            response: []

		}
	},
	methods: {
        // ricerca con "stringa input" and "type": taxomID or scientific name
		async nucleotideSearch() {
			this.loading = true;
			this.errormsg = null;
			try {
				//  
				let response = await this.$axios.get("/taxon_term", 
                {
                    params: 
                    {
                        search: this.search,
                        type: this.type
                    }
                    
				});
                this.response=response.data
			
			} catch (e) {
				this.errormsg = e.response.data
        
				}

			this.loading = false;
		},
	},


}

</script>

<template>
  <div>

    <div class="input">
      <label for="search-input">Input here</label>
      <input type="text" id="search-input" v-model="search" />
    </div>

    <div class="dropdown">
      <label for="dropdown">Type:</label>
      <select id="dropdown" v-model="type">
        <option value="scientific_name">Scientific Name</option>
        <option value="id">Taxom Id</option>
      </select>
    </div>

  </div>
    <button type="button" @click="nucleotideSearch()"> Search</button>


<!--TABELLA -->
    <ErrorMsg v-if="errormsg" :msg="errormsg"></ErrorMsg>
    <div>
    <table id="table">
      <thead>
        <tr>
          <th class="title" colspan="10">Datas</th>
        </tr>
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
            <a :href="'/organism/'+item.TaxId + '/taxon'">
              {{ item.ScientificName }}
            </a>
          </td>

          <td>
            <a :href="'/organism/'+item.TaxId+'/nucleotides'">

              {{ item.QtyNucleotides === 9999 ? item.QtyNucleotides+"+" : item.QtyNucleotides}}
            </a>
          </td>

          <td>
            <a :href="'/organism/'+item.TaxId + '/proteins'">
              {{ item.QtyProteins === 9999 ? item.QtyProteins+"+" : item.QtyProteins}}
            </a>
          </td>

          <td>
            <a :href="'/organism/'+item.TaxId + '/products'">
              {{ item.QtyProducts === 9999 ? item.QtyProducts+"+" : item.QtyProducts}}
            </a>
          </td>
            
          <td>{{ item.Genomes }}</td>
          <td>{{ item.Annotations }}</td>
          <td>{{ item.Trascriptome }}</td>
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
  padding: 10px 20px;
  background-color: #0d1aad;
  color: white;
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
   
#table a{
    color: #000;
}
   
</style>
