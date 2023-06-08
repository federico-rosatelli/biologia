<script>
import LoadingSpinner from '../components/LoadingSpinner.vue';


export default{
    data: function () {
        return {
            errormsg: null,
            loading: false,
            genomes: null,
            project: null,
        };
    },
    methods: {

        close(taxId){
            console.log(taxId);
            this.project = null
            document.getElementById("popup-"+taxId).style.display = "none"

        },

        ProjectsTree(genome) {
            if (this.genomeClick){
                this.close(TaxId)
            }
            this.genomeClick=genome
            document.getElementById("popup-"+genome.TaxId).style.display = "flex"
        },
        openModal(taxId){
          if (this.project){
            return
          }
          document.getElementById("popup-"+taxId).style.display = "flex"
          this.project = taxId
        },

        plus(count, add) {
          switch(count) {
            case ('p'):
              this.projcount += add
            case('s'):
              this.sampcount += add
            case('e'):
              this.expcount += add
            case('r'):
              this.runcount += add
          }

        }
    },
    async mounted(){
            try{
                let response = await this.$axios.get(`/analysis`);
                let projects = response.data
                let gen = {}
                projects.forEach(proj => {
                  if(!gen[proj.TaxId]){
                    gen[proj.TaxId] = []
                  }
                  gen[proj.TaxId].push(proj)
                });
                console.log(gen);
                this.genomes = gen
            } catch(e){
              console.log(e);
                this.errormsg = e.response
            }
        },
}


</script>

<template>
<br/>
<br/>
<br/>
<br/>
  <LoadingSpinner :loading="this.loading"/>
  <ErrorMsg v-if="errormsg" :msg="errormsg"></ErrorMsg>
  <div v-if="this.genomes">
    <div>
      <h1>{{ Object.keys(this.genomes).length }} Results</h1>
    </div>
    <table id="table">
      <thead>
        <tr>
          <th>Scientific Name</th>
          <th>Number BioProjects</th>
        </tr>
      </thead>
      <tr v-for="genome in this.genomes" :key="genome.TaxId">

        <td class="clickable">
          <div @click="openModal(genome[0].TaxId)">
            {{ genome[0].ScientificName ? genome[0].ScientificName : genome[0].TaxId }}
          </div>
          <div class="popup" :id="'popup-'+genome[0].TaxId">
            <div v-if="this.project" class="popup-content">
              <BioProjComp :data="genome"/>
              <button type="button" class="popupclose" @click="close(genome[0].TaxId)">Close</button>
            </div>
          </div>
        </td>
        <td>
          {{ genome.length }}
        </td>
        <!-- <div>
            <div class="popup" :id='"popup-"+genome.TaxId'>
                <div v-if="this.genomeClick" class="popup-content">
                    <BioProjComp :data="this.genomeClick"/>
                    <button type="button" class="popupclose" @click="close(genome.TaxId)">Close</button>
                </div>
            </div>
        </div>
              
          <td>
            {{ plus('p', genome.length) }}
            {{ projcount == 1 ? genome.BioProjects[0].BioProjectId : projcount }}
          </td>

          <td v-for="bioProject in genome" :key="bioProject.BioProjectId">
            {{ plus('s', bioProject.BioSamples.length) }}
            {{ sampcount == 1 ? bioProject.BioSamples[0].BioSampleId : sampcount }}
          </td>
          
          <td v-for="bioProject in genome" :key="bioProject.BioProjectId">
            <span v-for="bioSample in bioProject.BioSamples" :key="bioSample.BioSampleId">
              {{ plus('e', bioSample.Experiments.length) }}
              {{ expcount == 1 ? bioSample.Experiments[0].ExperimentId : expcount }}
            </span>
          </td>

          <td v-for="bioProject in genome" :key="bioProject.BioProjectId">
            <span v-for="bioSample in bioProject.BioSamples" :key="bioSample.BioSampleId">
              <span v-for="experiment in bioSample.Experiments" :key="experiment.ExperimentId">
                {{ plus('r', experiment.Runs.length) }}
                {{ runcount == 1 ? experiment.Runs[0] : runcount }}
              </span>
            </span>
          </td> -->

      </tr>
    </table>
  </div>
            


</template>

<style>
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

.popup{
    background: rgba(0, 0, 0, 0.6);
    position:absolute;
    /* overflow: auto;
    max-height: 100%; */
    top: 10%;
    left: 18%;
    right: 5%;
    /* right: 5%; */
    display: none;
    justify-content: center;
    align-items: center;
    border: 2px solid black;
}

.popup-content{
    overflow-y: auto;
    word-wrap: break-word;
    overflow-x:visible;
    top: 20%;
    bottom: 20%;
    right: 10%;
    left: 10%;
    width: 80%;
    color: #ffffff;
    size: 10px;
    box-shadow: 0 0 5px #0099FF;
    background:rgba(0, 0, 0, 0.815);
    padding: 2px;
    border-radius: 5px;
    position:fixed;
}


.popupclose{
  position: relative;
  bottom: 0px;
  right: 0px;
  border-radius: 5px;
  background-color: #003399;
}


.popup-close{
    position: absolute;
    top: -15px;
    right: -15px;
    background: #fff;
    height: 20px;
    width: 20px;
    border-radius: 50%;
    box-shadow: 6px 6px 29px -4px rgba(0,0,0,0.75);
    cursor: pointer;
}
</style>