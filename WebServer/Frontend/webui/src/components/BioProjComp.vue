<script>
export default {
	props: ['data'],
    methods:{
    }
}
</script>

<template>
    <!-- <div>BioProjects of {{ scientificName }}</div> -->
    <tr>
        BioProject
    </tr>
    <ul>
        <tr v-for="proj in data">
            BioProject: {{ proj.BioProjectId }}
            <td v-for="sample in proj.BioSamples">
                BioSample: {{ sample.BioSampleId }}
                <td v-for="exp in sample.Experiments">
                    Experiment: {{ exp.ExperimentId }}
                    <ul v-for="run in exp.Runs">
                        Run: {{ run }}
                    </ul>
                </td>
            </td>
        </tr>
    </ul>
    <div v-if="data[0].ImagesPath && data[0].ImagesPath.length > 0" class="immagini">
        <div v-for="img in data[0].ImagesPath">
            <img :src="img" alt="Prova">
        </div>
    </div>
</template>