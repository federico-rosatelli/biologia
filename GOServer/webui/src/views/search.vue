<!--IDEA: Vue delle ricerche con funzine che prende tipo di ricerca e input e ritorna tabella -->
<script>

export default {
	data: function() {
		return {
			errormsg: null,
			loading: false,
      search: null,
      type: "scientific_name",
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
    <div class="dropdown">
      <label for="dropdown">Type:</label>
      <select id="dropdown" name="drop" v-model="type">
        <option value="scientific_name" selected>Scientific Name</option>
        <option value="id">Taxom Id</option>
      </select>
    </div>
    <div class="input">
      <input type="text" id="search-input" placeholder="Input here" v-model="search" />
      <button type="button" @click="nucleotideSearch()"> Search</button>
    </div>
  </div>


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
            <RouterLink :to="'/organism/'+item.TaxId + '/taxon'">
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
            <RouterLink :to="'/organism/'+item.TaxId + '/products'">
              {{ item.QtyProducts >= 9999 ? 9999+"+" : item.QtyProducts}}
            </RouterLink>
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
   
</style>