<script>
import LoadingSpinner from '../components/LoadingSpinner.vue'

export default {
    data: function () {
        return {
            errormsg: null,
            loading: false,
            runData: $axios.params.runDataId,
            runDataType: null
        };
    },
    methods: {

        async mounted(){
            try{
                let response = await this.$axios.get(`/sra` + this.runData);
                this.runDataType=response.data
                
            } catch(e){
                this.errormsg = e.response.data
            }
        },
    },
    components: { LoadingSpinner }
}

</script>

<template>
    <div>
      <h2>Images:</h2>
      <ul>
        <li v-for="imagePath in runDataType.ImagesPath" :key="imagePath">
          <img :src="imagePath" alt="Image">
        </li>
      </ul>
  
      <h2>CSV File:</h2>
      <ul>
        <li v-for="(values, fieldName) in runDataType.CsvFile" :key="fieldName">
          <strong>{{ fieldName }}:</strong>
          <ul>
            <li v-for="value in values" :key="value">
                {{ value }}
            </li>
          </ul>
        </li>
      </ul>
    </div>
  </template>
  