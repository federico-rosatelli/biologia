<!--IDEA: Vue delle ricerche con funzine che prende tipo di ricerca e input e ritorna tabella -->
<script>

export default {
	data: function() {
		return {
			errormsg: null,
			loading: false,
            search: null,
            type: null,
            response: [],
            nucleo:null,
            dataArray: [],
            otherRouterLink:null,
            taxonomy:null

		}
	},
    
	methods: {
        async findLocus(locus){
            this.loading = true
            if (this.nucleo){
                this.close(this.nucleo.GBSeqLocus)
            }
            try{
                let taxId = this.$route.params.taxid;
                let typof = this.$route.path.split("/").pop();
                let nucleo = await this.$axios.get(`organism/${taxId}/${typof}/${locus}`);
                this.nucleo= typof == "proteins" ? nucleo.data.Proteins[0] : nucleo.data.Nucleotides[0]
            } catch(e){
                this.errormsg = e.nucleo
            }
            document.getElementById("popup-"+locus).style.display = "flex"
            this.loading = false
        },
        close(locus){
            document.getElementById("popup-"+locus).style.display = "none"
            this.nucleo = null
        }
    },
    async mounted(){
        let taxId = this.$route.params.taxid;
        let typof = this.$route.path.split("/").pop();
        try{
            let response = await this.$axios.get(`/organism/${taxId}/${typof}`);
            this.response=response.data
        } catch(e){
            this.errormsg = e.response.data
        }
        if (!this.errormsg){
            try{
                let taxonomy = await this.$axios.get(`/taxonomy?search=${this.response.ScientificName}`);
                this.taxonomy=taxonomy.data
            } catch(e){
                this.errormsg = e.response.data
            }
        }

        this.dataArray = typof == "proteins" ? this.response.Proteins : this.response.Nucleotides
        let datalink = this.$route.path.split('/').pop() == 'nucleotides' ? '/proteins': '/nucleotides'
        this.otherRouterLink = '/organism/'+this.$route.params.taxid+datalink

    },
}


</script>

<template>
<!--TABELLA -->
    <div>
        <a v-if="otherRouterLink" :href="otherRouterLink">{{otherRouterLink}}</a>
    </div>
    <ErrorMsg v-if="errormsg" :msg="errormsg"></ErrorMsg>
    <div v-if="this.dataArray">
        <br/>
        <LoadingSpinner :loading="this.loading"/>
        <div style="display: flex;">

            <div class="title-left">
                Locus Id
            </div>
            <div class="title-right">
                Taxonomy
            </div>
        </div>
        <div class="block-organism">
            <div class="left">
                <ul>
                    <div v-for="item in this.dataArray">
                        <div @click="findLocus(item.GBSeq_locus)" class="clickable">
                            {{ item.GBSeq_locus }}
                        </div>
                        <div>
                            <div class="popup" :id='"popup-"+item.GBSeq_locus'>
                                <div v-if="this.nucleo" class="popup-content">
                                    <NucleoComp :data="this.nucleo"/>
                                    <button type="button" class="popupclose" @click="close(item.GBSeq_locus)">Close</button>
                                </div>
                            </div>
                            
                        </div>
                    </div>
                </ul>
            </div>
            <div class="right" v-if="taxonomy">
                <h2>
                    Organism Name: {{ taxonomy.ScientificName }}
                </h2>
                <h2>
                    Taxonomy Id: {{ taxonomy.TaxID }}
                </h2>
                <h2>
                    Rank: {{ taxonomy.Rank }}
                </h2>
                <h2>
                    Division: {{ taxonomy.Division }}
                </h2>
                <h2>
                    Lineage:
                    <ul v-for="lin in taxonomy.LineageEx">
                        <h3>
                            -----------------------------
                        </h3>
                        <a :href="'/taxonomy/'+lin.TaxID">
                            <h3>
                                Taxonomy Id: {{ lin.TaxID }}
                            </h3>
                            <h3>
                                Scientific Name: {{ lin.ScientificName }}
                            </h3>
                            <h3>
                                Rank: {{ lin.Rank }}
                            </h3>
                        </a>
                        <h3>
                            -----------------------------
                        </h3>
                    </ul>
                </h2>
                
            </div>
        </div>
    </div>
    

</template>

<style>
.clickable:hover, .clickable:focus, .clickable:active { color:#ffcc00; }

.block-organism {
    display: flex;
}
.title-left{
    margin-right: auto;
    font-weight: bold;
}
.title-right{
    margin-left: auto;
    font-weight: bold;
    margin-right: 33%;
}

.block-organism .left{
    background: rgba(0, 0, 0, 0.6);
    /* overflow: auto;
    max-height: 100%; */
    width: 20%;
    max-height: 500px;
    overflow-y: auto;
    justify-content: center;
    align-items: center;
    border: 2px solid black;
}
.block-organism .right {
    width: 40%;
    margin-left: auto;
    background: rgba(0, 0, 0, 0.6);
    /* overflow: auto;
    max-height: 100%; */
    justify-content: center;
    align-items: center;
    border: 2px solid black;
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
