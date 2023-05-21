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
            nucleo:null

		}
	},
    
	methods: {
        async findLocus(locus){
            if (this.nucleo){
                this.close(this.nucleo.GBSeqLocus)
            }
            try{
                let taxId = this.$route.params.taxid;
                let nucleo = await this.$axios.get(`organism/${taxId}/nucleotide/${locus}`);
                this.nucleo=nucleo.data.Nucleotides[0]
            } catch(e){
                this.errormsg = e.response.data
            }
            document.getElementById("popup-"+locus).style.display = "flex"
        },
        close(locus){
            document.getElementById("popup-"+locus).style.display = "none"
            this.nucleo = null
        }
    },
    async mounted(){
        let taxId = this.$route.params.taxid;
        try{
            let response = await this.$axios.get(`/organism/${taxId}/nucleotides`);
            this.response=response.data
        } catch(e){
            this.errormsg = e.response.data
        }
    },
}


</script>

<template>
<!--TABELLA -->
   
    <ErrorMsg v-if="errormsg" :msg="errormsg"></ErrorMsg>
    <div v-else>
        <div class="clickable">
            <h2>
                {{ response.ScientificName }}
            </h2>
        </div>
        <ul>
            <div v-for="item in response.Nucleotides">
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


</template>

<style>
.clickable:hover, .clickable:focus, .clickable:active { color:#ffcc00; }
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
    overflow-y: scroll;
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
