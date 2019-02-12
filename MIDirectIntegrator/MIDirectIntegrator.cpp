#include <mitsuba/render/scene.h>

MTS_NAMESPACE_BEGIN
/*
used to test important sample envmap and bsdf!
*/
class MIDirectIntegrator : public SamplingIntegrator 
{
public:
	Float m_bsdf_frac;
	bool m_sample_sum;
	bool m_change_bsdf_fraction;
	bool m_change_bsdf_fraction_use_fix_fraction;
	MIDirectIntegrator(const Properties &props) : SamplingIntegrator(props) 
	{
		/* Number of shading samples -- this parameter is a shorthand notation
		to set both 'emitterSamples' and 'bsdfSamples' at the same time*/
		size_t shadingSamples = props.getSize("shadingSamples", 1);

		/* Number of samples to take using the emitter sampling technique */
		m_emitterSamples = props.getSize("emitterSamples", shadingSamples);
		/* Number of samples to take using the BSDF sampling technique */
		m_bsdfSamples = props.getSize("bsdfSamples", shadingSamples);
		/* Be strict about potential inconsistencies involving shading normals? */
		m_strictNormals = props.getBoolean("strictNormals", false);
		/* When this flag is set to true, contributions from directly
		* visible emitters will not be included in the rendered image */
		m_hideEmitters = props.getBoolean("hideEmitters", false);
		m_bsdf_frac = props.getFloat("bsdfFraction",0.5f);
		m_sample_sum = props.getBoolean("sampleSum",false);
		m_change_bsdf_fraction = props.getBoolean("changeBsdfFracthion", false);
		m_change_bsdf_fraction_use_fix_fraction = props.getBoolean("changeBsdfFracthion_useFixFraction", false);
		
		Assert(m_emitterSamples + m_bsdfSamples > 0);
	}

	/// Unserialize from a binary data stream
	MIDirectIntegrator(Stream *stream, InstanceManager *manager)
		: SamplingIntegrator(stream, manager) 
	{
		m_emitterSamples = stream->readSize();
		m_bsdfSamples = stream->readSize();
		m_strictNormals = stream->readBool();
		m_hideEmitters = stream->readBool();
		configure();
	}

	void serialize(Stream *stream, InstanceManager *manager) const 
	{
		SamplingIntegrator::serialize(stream, manager);
		stream->writeSize(m_emitterSamples);
		stream->writeSize(m_bsdfSamples);
		stream->writeBool(m_strictNormals);
		stream->writeBool(m_hideEmitters);
	}

	void configure() 
	{
		SamplingIntegrator::configure();

		size_t sum = m_emitterSamples + m_bsdfSamples;
		m_weightBSDF = 1 / (Float)m_bsdfSamples;
		m_weightLum = 1 / (Float)m_emitterSamples;
		m_fracBSDF = m_bsdfSamples / (Float)sum;
		m_fracLum = m_emitterSamples / (Float)sum;
	}

	void configureSampler(const Scene *scene, Sampler *sampler) 
	{
		SamplingIntegrator::configureSampler(scene, sampler);
		if (m_emitterSamples > 1)
			sampler->request2DArray(m_emitterSamples);
		if (m_bsdfSamples > 1)
			sampler->request2DArray(m_bsdfSamples);
	}

	Spectrum LiChange(const RayDifferential &r, RadianceQueryRecord &rRec) const
	{
		/* Some aliases and local variables */
		const Scene *scene = rRec.scene;
		Intersection &its = rRec.its;
		RayDifferential ray(r);
		Spectrum Li(0.0f);
		Point2 sample;

		/* Perform the first ray intersection (or ignore if the
		intersection has already been provided). */
		if (!rRec.rayIntersect(ray)) {
			/* If no intersection could be found, possibly return
			radiance from a background emitter */
			if (rRec.type & RadianceQueryRecord::EEmittedRadiance && !m_hideEmitters)
				return scene->evalEnvironment(ray);
			else
				return Spectrum(0.0f);
		}

		/* Possibly include emitted radiance if requested */
		if (its.isEmitter() && (rRec.type & RadianceQueryRecord::EEmittedRadiance) && !m_hideEmitters)
			Li += its.Le(-ray.d);

		/* Include radiance from a subsurface scattering model if requested */
		if (its.hasSubsurface() && (rRec.type & RadianceQueryRecord::ESubsurfaceRadiance))
			Li += its.LoSub(scene, rRec.sampler, -ray.d, rRec.depth);

		const BSDF *bsdf = its.getBSDF(ray);

		if (!(rRec.type & RadianceQueryRecord::EDirectSurfaceRadiance)
			|| (m_strictNormals && dot(ray.d, its.geoFrame.n)
				* Frame::cosTheta(its.wi) >= 0)) {
			/* Only render the direct illumination component if
			*
			* 1. It was requested
			* 2. The surface has an associated BSDF (i.e. it isn't an index-
			*    matched medium transition -- this is not supported by 'direct')
			* 3. If 'strictNormals'=true, when the geometric and shading
			*    normals classify the incident direction to the same side
			*/
			return Li;
		}

		/* ==================================================================== */
		/*                          Emitter sampling                          */
		/* ==================================================================== */
		bool adaptiveQuery = (rRec.extra & RadianceQueryRecord::EAdaptiveQuery);

		/* Figure out how many BSDF and direct illumination samples to
		generate, and where the random numbers should come from */
		Point2 *sampleArray;
		size_t numDirectSamples = m_emitterSamples,
			numBSDFSamples = m_bsdfSamples;
		Float fracLum = m_fracLum, fracBSDF = m_fracBSDF,
			weightLum = m_weightLum, weightBSDF = m_weightBSDF;

		if (rRec.depth > 1 || adaptiveQuery) {
			/* This integrator is used recursively by another integrator.
			Be less accurate as this sample will not directly be observed. */
			numBSDFSamples = numDirectSamples = 1;
			fracLum = fracBSDF = .5f;
			weightLum = weightBSDF = 1.0f;
		}

		if (numBSDFSamples > 1) {
			sampleArray = rRec.sampler->next2DArray(numBSDFSamples);
		}
		else {
			sample = rRec.nextSample2D(); sampleArray = &sample;
		}

		Spectrum li_sample_light(0.f), li_sample_bsdf(0.f);
		Spectrum test_li_sample_light(0.f), test_li_sample_bsdf(0.f);
		Float pdf_sample_light=0.f, pdf_sample_bsdf=0.f;
		Float pdf_eval_light=0.f, pdf_eval_bsdf=0.f;

		DirectSamplingRecord dRec(its);
		if (bsdf->getType() & BSDF::ESmooth)
		{
			/* Only use direct illumination sampling when the surface's
			BSDF has smooth (i.e. non-Dirac delta) component */
			for (size_t i = 0; i<1; ++i) {
				/* Estimate the direct illumination if this is requested */
				Spectrum value = scene->sampleEmitterDirect(dRec, sampleArray[i]);
				pdf_sample_light = dRec.pdf;
				if (!value.isZero()) 
				{
					
					value *= pdf_sample_light;
					const Emitter *emitter = static_cast<const Emitter *>(dRec.object);

					/* Allocate a record for querying the BSDF */
					BSDFSamplingRecord bRec(its, its.toLocal(dRec.d));

					/* Evaluate BSDF * cos(theta) */
					const Spectrum bsdfVal = bsdf->eval(bRec);

					if (!bsdfVal.isZero() && (!m_strictNormals
						|| dot(its.geoFrame.n, dRec.d) * Frame::cosTheta(bRec.wo) > 0)) {
						/* Calculate prob. of sampling that direction using BSDF sampling */
						Float bsdfPdf = emitter->isOnSurface() ? bsdf->pdf(bRec) : 0;
						pdf_eval_bsdf = bsdfPdf;
						test_li_sample_light = value / pdf_sample_light*bsdfVal;
						li_sample_light = value*bsdfVal;
					}
					else
						pdf_eval_bsdf = 0.f;
				}
				else
				{
					pdf_eval_bsdf = 0.f;
				}
				//else return Li;
			}
		}

		/* ==================================================================== */
		/*                            BSDF sampling                             */
		/* ==================================================================== */

		if (numBSDFSamples > 1) {
			sampleArray = rRec.sampler->next2DArray(numBSDFSamples);
		}
		else {
			sample = rRec.nextSample2D(); sampleArray = &sample;
		}

		Intersection bsdfIts;
		for (size_t i = 0; i<1; ++i) {
			/* Sample BSDF * cos(theta) and also request the local density */
			Float bsdfPdf;

			BSDFSamplingRecord bRec(its, rRec.sampler, ERadiance);
			Spectrum bsdfVal = bsdf->sample(bRec, bsdfPdf, sampleArray[i]);
			if (bsdfVal.isZero())
				break;// return Li;
			pdf_sample_bsdf = bsdfPdf;
			bsdfVal *= bsdfPdf;
			/* Prevent light leaks due to the use of shading normals */
			const Vector wo = its.toWorld(bRec.wo);
			Float woDotGeoN = dot(its.geoFrame.n, wo);
			if (m_strictNormals && woDotGeoN * Frame::cosTheta(bRec.wo) <= 0)
				break;// return Li;

			/* Trace a ray in this direction */
			Ray bsdfRay(its.p, wo, ray.time);

			Spectrum value;
			if (scene->rayIntersect(bsdfRay, bsdfIts)) {
				/* Intersected something - check if it was an emitter */
				if (!bsdfIts.isEmitter())
					break;// return Li;

				value = bsdfIts.Le(-bsdfRay.d);
				dRec.setQuery(bsdfRay, bsdfIts);
			}
			else {
				/* Intersected nothing -- perhaps there is an environment map? */
				const Emitter *env = scene->getEnvironmentEmitter();

				if (!env || (m_hideEmitters && bRec.sampledType == BSDF::ENull))
					break;// return Li;

				value = env->evalEnvironment(RayDifferential(bsdfRay));
				if (!env->fillDirectSamplingRecord(dRec, bsdfRay))
					break;// return Li;
			}

			/* Compute the prob. of generating that direction using the
			implemented direct illumination sampling technique */
			const Float lumPdf = (!(bRec.sampledType & BSDF::EDelta)) ?
				scene->pdfEmitterDirect(dRec) : 0;
			pdf_eval_light = lumPdf;
			test_li_sample_bsdf = bsdfVal / bsdfPdf*value;
			li_sample_bsdf = value*bsdfVal;
		}
		if ((test_li_sample_bsdf.average() + test_li_sample_light.average()) == 0.f)
			return Li;
		Float frc = m_change_bsdf_fraction_use_fix_fraction?m_bsdf_frac:
			test_li_sample_bsdf.average() / (test_li_sample_bsdf.average() + test_li_sample_light.average());
		//if (pdf_sample_light <= 0.f || pdf_eval_bsdf == 0.f||
		if(rRec.nextSample1D() < frc)
		{
			if ((frc*pdf_sample_bsdf + (1 - frc)*pdf_eval_light) == 0.f)
				return Li;
			Li += li_sample_bsdf / (frc*pdf_sample_bsdf+(1-frc)*pdf_eval_light);
		}
		else
		{
			if ((frc*pdf_eval_bsdf + (1 - frc)*pdf_sample_light) == 0.f)
				return Li;
			Li += li_sample_light / (frc*pdf_eval_bsdf + (1 - frc)*pdf_sample_light);
		}

		return Li;
		
#if 0
		if (bsdf->getType() & BSDF::ESmooth) {
			/* Only use direct illumination sampling when the surface's
			BSDF has smooth (i.e. non-Dirac delta) component */
			for (size_t i = 0; i<numDirectSamples; ++i) {
				/* Estimate the direct illumination if this is requested */
				Spectrum value = scene->sampleEmitterDirect(dRec, sampleArray[i]);
				if (!value.isZero()) {
					const Emitter *emitter = static_cast<const Emitter *>(dRec.object);

					/* Allocate a record for querying the BSDF */
					BSDFSamplingRecord bRec(its, its.toLocal(dRec.d));

					/* Evaluate BSDF * cos(theta) */
					const Spectrum bsdfVal = bsdf->eval(bRec);

					if (!bsdfVal.isZero() && (!m_strictNormals
						|| dot(its.geoFrame.n, dRec.d) * Frame::cosTheta(bRec.wo) > 0)) {
						/* Calculate prob. of sampling that direction using BSDF sampling */
						Float bsdfPdf = emitter->isOnSurface() ? bsdf->pdf(bRec) : 0;

						/* Weight using the power heuristic */
						const Float weight = miWeight(dRec.pdf * fracLum,
							bsdfPdf * fracBSDF) * weightLum;

						Li += value * bsdfVal * weight;
					}
				}
			}
		}

		/* ==================================================================== */
		/*                            BSDF sampling                             */
		/* ==================================================================== */

		if (numBSDFSamples > 1) {
			sampleArray = rRec.sampler->next2DArray(numBSDFSamples);
		}
		else {
			sample = rRec.nextSample2D(); sampleArray = &sample;
		}

		Intersection bsdfIts;
		for (size_t i = 0; i<numBSDFSamples; ++i) {
			/* Sample BSDF * cos(theta) and also request the local density */
			Float bsdfPdf;

			BSDFSamplingRecord bRec(its, rRec.sampler, ERadiance);
			Spectrum bsdfVal = bsdf->sample(bRec, bsdfPdf, sampleArray[i]);
			if (bsdfVal.isZero())
				continue;

			/* Prevent light leaks due to the use of shading normals */
			const Vector wo = its.toWorld(bRec.wo);
			Float woDotGeoN = dot(its.geoFrame.n, wo);
			if (m_strictNormals && woDotGeoN * Frame::cosTheta(bRec.wo) <= 0)
				continue;

			/* Trace a ray in this direction */
			Ray bsdfRay(its.p, wo, ray.time);

			Spectrum value;
			if (scene->rayIntersect(bsdfRay, bsdfIts)) {
				/* Intersected something - check if it was an emitter */
				if (!bsdfIts.isEmitter())
					continue;

				value = bsdfIts.Le(-bsdfRay.d);
				dRec.setQuery(bsdfRay, bsdfIts);
			}
			else {
				/* Intersected nothing -- perhaps there is an environment map? */
				const Emitter *env = scene->getEnvironmentEmitter();

				if (!env || (m_hideEmitters && bRec.sampledType == BSDF::ENull))
					continue;

				value = env->evalEnvironment(RayDifferential(bsdfRay));
				if (!env->fillDirectSamplingRecord(dRec, bsdfRay))
					continue;
			}

			/* Compute the prob. of generating that direction using the
			implemented direct illumination sampling technique */
			const Float lumPdf = (!(bRec.sampledType & BSDF::EDelta)) ?
				scene->pdfEmitterDirect(dRec) : 0;

			/* Weight using the power heuristic */
			const Float weight = miWeight(bsdfPdf * fracBSDF,
				lumPdf * fracLum) * weightBSDF;

			Li += value * bsdfVal * weight;
		}

		return Li;
#endif
	}

	Spectrum Li(const RayDifferential &r, RadianceQueryRecord &rRec) const 
	{
		if (m_change_bsdf_fraction)
		{
			return LiChange(r,rRec);
		}
		/* Some aliases and local variables */
		const Scene *scene = rRec.scene;
		Intersection &its = rRec.its;
		RayDifferential ray(r);
		Spectrum Li(0.0f);
		Point2 sample;

		/* Perform the first ray intersection (or ignore if the
		intersection has already been provided). */
		if (!rRec.rayIntersect(ray)) 
		{
			/* If no intersection could be found, possibly return
			radiance from a background emitter */
			if (rRec.type & RadianceQueryRecord::EEmittedRadiance && !m_hideEmitters)
				return scene->evalEnvironment(ray);
			else
				return Spectrum(0.0f);
		}

		/* Possibly include emitted radiance if requested */
		if (its.isEmitter() && (rRec.type & RadianceQueryRecord::EEmittedRadiance) && !m_hideEmitters)
			Li += its.Le(-ray.d);

		/* Include radiance from a subsurface scattering model if requested */
		if (its.hasSubsurface() && (rRec.type & RadianceQueryRecord::ESubsurfaceRadiance))
			Li += its.LoSub(scene, rRec.sampler, -ray.d, rRec.depth);

		const BSDF *bsdf = its.getBSDF(ray);

		if (!(rRec.type & RadianceQueryRecord::EDirectSurfaceRadiance)
			|| (m_strictNormals && dot(ray.d, its.geoFrame.n)
				* Frame::cosTheta(its.wi) >= 0)) 
		{
			/* Only render the direct illumination component if
			*
			* 1. It was requested
			* 2. The surface has an associated BSDF (i.e. it isn't an index-
			*    matched medium transition -- this is not supported by 'direct')
			* 3. If 'strictNormals'=true, when the geometric and shading
			*    normals classify the incident direction to the same side
			*/
			return Li;
		}

		/* ==================================================================== */
		/*                          Emitter sampling                          */
		/* ==================================================================== */
		bool adaptiveQuery = (rRec.extra & RadianceQueryRecord::EAdaptiveQuery);

		/* Figure out how many BSDF and direct illumination samples to
		generate, and where the random numbers should come from */
		Point2 *sampleArray;
		size_t numDirectSamples = m_emitterSamples,
			numBSDFSamples = m_bsdfSamples;
		Float fracLum = m_fracLum, fracBSDF = m_fracBSDF,
			weightLum = m_weightLum, weightBSDF = m_weightBSDF;

		if (rRec.depth > 1 || adaptiveQuery) 
		{
			/* This integrator is used recursively by another integrator.
			Be less accurate as this sample will not directly be observed. */
			numBSDFSamples = numDirectSamples = 1;
			fracLum = fracBSDF = .5f;
			weightLum = weightBSDF = 1.0f;
		}

		if (numDirectSamples > 1) 
		{
			sampleArray = rRec.sampler->next2DArray(numDirectSamples);
		}
		else 
		{
			sample = rRec.nextSample2D(); sampleArray = &sample;
		}

		if (!m_sample_sum)
		{
			DirectSamplingRecord dRec(its);
			if (bsdf->getType() & BSDF::ESmooth) 
			{
				/* Only use direct illumination sampling when the surface's
				BSDF has smooth (i.e. non-Dirac delta) component */
				for (size_t i = 0; i<numDirectSamples; ++i) {
					/* Estimate the direct illumination if this is requested */
					Spectrum value = scene->sampleEmitterDirect(dRec, sampleArray[i]);
					if (!value.isZero()) {
						const Emitter *emitter = static_cast<const Emitter *>(dRec.object);

						/* Allocate a record for querying the BSDF */
						BSDFSamplingRecord bRec(its, its.toLocal(dRec.d));

						/* Evaluate BSDF * cos(theta) */
						const Spectrum bsdfVal = bsdf->eval(bRec);

						if (!bsdfVal.isZero() && (!m_strictNormals
							|| dot(its.geoFrame.n, dRec.d) * Frame::cosTheta(bRec.wo) > 0)) {
							/* Calculate prob. of sampling that direction using BSDF sampling */
							Float bsdfPdf = emitter->isOnSurface() ? bsdf->pdf(bRec) : 0;

							/* Weight using the power heuristic */
							const Float weight = miWeight(dRec.pdf * fracLum,bsdfPdf * fracBSDF) * weightLum;

							Li += value * bsdfVal * weight;
						}
					}
				}
			}
			
			/* ==================================================================== */
			/*                            BSDF sampling                             */
			/* ==================================================================== */

			if (numBSDFSamples > 1) {
				sampleArray = rRec.sampler->next2DArray(numBSDFSamples);
			}
			else {
				sample = rRec.nextSample2D(); sampleArray = &sample;
			}

			Intersection bsdfIts;
			for (size_t i = 0; i<numBSDFSamples; ++i) {
				/* Sample BSDF * cos(theta) and also request the local density */
				Float bsdfPdf;

				BSDFSamplingRecord bRec(its, rRec.sampler, ERadiance);
				Spectrum bsdfVal = bsdf->sample(bRec, bsdfPdf, sampleArray[i]);
				if (bsdfVal.isZero())
					continue;

				/* Prevent light leaks due to the use of shading normals */
				const Vector wo = its.toWorld(bRec.wo);
				Float woDotGeoN = dot(its.geoFrame.n, wo);
				if (m_strictNormals && woDotGeoN * Frame::cosTheta(bRec.wo) <= 0)
					continue;

				/* Trace a ray in this direction */
				Ray bsdfRay(its.p, wo, ray.time);

				Spectrum value;
				if (scene->rayIntersect(bsdfRay, bsdfIts)) {
					/* Intersected something - check if it was an emitter */
					if (!bsdfIts.isEmitter())
						continue;

					value = bsdfIts.Le(-bsdfRay.d);
					dRec.setQuery(bsdfRay, bsdfIts);
				}
				else {
					/* Intersected nothing -- perhaps there is an environment map? */
					const Emitter *env = scene->getEnvironmentEmitter();

					if (!env || (m_hideEmitters && bRec.sampledType == BSDF::ENull))
						continue;

					value = env->evalEnvironment(RayDifferential(bsdfRay));
					if (!env->fillDirectSamplingRecord(dRec, bsdfRay))
						continue;
				}

				/* Compute the prob. of generating that direction using the
				implemented direct illumination sampling technique */
				const Float lumPdf = (!(bRec.sampledType & BSDF::EDelta)) ?
					scene->pdfEmitterDirect(dRec) : 0;

				/* Weight using the power heuristic */
				const Float weight = miWeight(bsdfPdf * fracBSDF,lumPdf * fracLum) * weightBSDF;

				Li += value * bsdfVal * weight;
			}

			return Li;
		}
		else
		{

			Float sample_fraction = m_bsdf_frac;
			DirectSamplingRecord dRec(its);
			if (rRec.nextSample1D() >= sample_fraction)
			{


				if (bsdf->getType() & BSDF::ESmooth)
				{
					/* Only use direct illumination sampling when the surface's
					BSDF has smooth (i.e. non-Dirac delta) component */
					for (size_t i = 0; i < numDirectSamples; ++i)
					{
						/* Estimate the direct illumination if this is requested */
						Spectrum value = scene->sampleEmitterDirect(dRec, sampleArray[i]);
						if (!value.isZero())
						{
							const Emitter *emitter = static_cast<const Emitter *>(dRec.object);
							/* Allocate a record for querying the BSDF */
							BSDFSamplingRecord bRec(its, its.toLocal(dRec.d));
							Float pdf_emitter = dRec.pdf;
							value *= pdf_emitter;
							/* Evaluate BSDF * cos(theta) */
							const Spectrum bsdfVal = bsdf->eval(bRec);

							if (!bsdfVal.isZero() && (!m_strictNormals
								|| dot(its.geoFrame.n, dRec.d) * Frame::cosTheta(bRec.wo) > 0))
							{
								/* Calculate prob. of sampling that direction using BSDF sampling */
								Float bsdfPdf = emitter->isOnSurface() ? bsdf->pdf(bRec) : 0;

								/* Weight using the power heuristic */
								Float weight = miWeight(dRec.pdf * fracLum, bsdfPdf * fracBSDF) * weightLum;
								weight = 1.f / (bsdfPdf*sample_fraction + (1.f - sample_fraction)*pdf_emitter);

								Li += value * bsdfVal * weight;
							}
						}
					}
				}
			}
			else
			{
				/* ==================================================================== */
				/*                            BSDF sampling                             */
				/* ==================================================================== */

				if (numBSDFSamples > 1) {
					sampleArray = rRec.sampler->next2DArray(numBSDFSamples);
				}
				else {
					sample = rRec.nextSample2D(); sampleArray = &sample;
				}

				Intersection bsdfIts;
				for (size_t i = 0; i < numBSDFSamples; ++i) {
					/* Sample BSDF * cos(theta) and also request the local density */
					Float bsdfPdf;

					BSDFSamplingRecord bRec(its, rRec.sampler, ERadiance);
					Spectrum bsdfVal = bsdf->sample(bRec, bsdfPdf, sampleArray[i]);
					if (bsdfVal.isZero())
						continue;
					bsdfVal *= bsdfPdf;

					/* Prevent light leaks due to the use of shading normals */
					const Vector wo = its.toWorld(bRec.wo);
					Float woDotGeoN = dot(its.geoFrame.n, wo);
					if (m_strictNormals && woDotGeoN * Frame::cosTheta(bRec.wo) <= 0)
						continue;

					/* Trace a ray in this direction */
					Ray bsdfRay(its.p, wo, ray.time);

					Spectrum value;
					if (scene->rayIntersect(bsdfRay, bsdfIts))
					{
						/* Intersected something - check if it was an emitter */
						if (!bsdfIts.isEmitter())
							continue;

						value = bsdfIts.Le(-bsdfRay.d);
						dRec.setQuery(bsdfRay, bsdfIts);
					}
					else
					{
						/* Intersected nothing -- perhaps there is an environment map? */
						const Emitter *env = scene->getEnvironmentEmitter();

						if (!env || (m_hideEmitters && bRec.sampledType == BSDF::ENull))
							continue;

						value = env->evalEnvironment(RayDifferential(bsdfRay));
						if (!env->fillDirectSamplingRecord(dRec, bsdfRay))
							continue;
					}

					/* Compute the prob. of generating that direction using the
					implemented direct illumination sampling technique */
					const Float lumPdf = (!(bRec.sampledType & BSDF::EDelta)) ?
						scene->pdfEmitterDirect(dRec) : 0;

					/* Weight using the power heuristic */
					Float weight = miWeight(bsdfPdf * fracBSDF, lumPdf * fracLum) * weightBSDF;
					weight = 1.f / (bsdfPdf*sample_fraction + (1.f - sample_fraction)*lumPdf);

					Li += value * bsdfVal * weight;
				}
			}
			return Li;
		}
	}

	inline Float miWeight(Float pdfA, Float pdfB) const {
		pdfA *= pdfA; pdfB *= pdfB;
		return pdfA / (pdfA + pdfB);
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "MIDirectIntegrator[" << endl
			<< "  emitterSamples = " << m_emitterSamples << "," << endl
			<< "  bsdfSamples = " << m_bsdfSamples << "," << endl
			<< "  strictNormals = " << m_strictNormals << endl
			<< "]";
		return oss.str();
	}

	MTS_DECLARE_CLASS()
private:
	size_t m_emitterSamples;
	size_t m_bsdfSamples;
	Float m_fracBSDF, m_fracLum;
	Float m_weightBSDF, m_weightLum;
	bool m_strictNormals;
	bool m_hideEmitters;
};
#endif

MTS_IMPLEMENT_CLASS_S(MIDirectIntegrator, false, SamplingIntegrator)
MTS_EXPORT_PLUGIN(MIDirectIntegrator, "Direct illumination integrator");
MTS_NAMESPACE_END
