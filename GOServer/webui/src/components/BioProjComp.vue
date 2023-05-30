<script>
export default {
	props: ['data'],
    methods:{

    }
}
</script>

<template>
    <div>BioProjects of {{ data.ScientificName }}</div>

    <ul v-for="proj in data.BioProjects" :key=proj.BioProjectId>
        {{ proj.BioProjectId }} project

        <ul v-for="bioProject in data.BioProjects" :key="bioProject.BioProjectId">
            <span v-for="bioSample in bioProject.BioSamples" :key="bioSample.BioSampleId">
                {{ bioSample.BioSampleId }} sample
            </span>
        </ul>

        <ul v-for="bioProject in data.BioProjects" :key="bioProject.BioProjectId">
            <span v-for="bioSample in bioProject.BioSamples" :key="bioSample.BioSampleId">
                <span v-for="experiment in bioSample.Experiments" :key="experiment.ExperimentId">
                    {{ experiment.ExperimentId }} experiment
                </span>
            </span>
        </ul>

        <ul v-for="bioProject in data.BioProjects" :key="bioProject.BioProjectId">
            <span v-for="bioSample in bioProject.BioSamples" :key="bioSample.BioSampleId">
                <span v-for="experiment in bioSample.Experiments" :key="experiment.ExperimentId">
                    <span v-for="run in experiment.Runs" :key="run.RunId">
                        {{ run.RunId }} run

                    </span>
                </span>
            </span>
        </ul>

        <ul v-for="bioProject in data.BioProjects" :key="bioProject.BioProjectId">
            <span v-for="bioSample in bioProject.BioSamples" :key="bioSample.BioSampleId">
                <span v-for="experiment in bioSample.Experiments" :key="experiment.ExperimentId">
                    <span v-for="run in experiment.Runs" :key="run.RunId">
                        <span v-for="rundata in run.RunDataId" :key="rundata.RunDataIdValue">
                            concorrenza con {{ rundata.ConcorrenzaCon }}
                            <a :href="'/sra/'+ rundata.RunDataIdValue"> <!-- non so quale sarà l'url preciso-->
                                {{ rundata.RunDataIdValue }}
                            </a>

                        </span>

                    </span>
                </span>
            </span>
        </ul>

    </ul>
</template>