

CREATE TABLE "d_patients" (
	"subject_id" integer NOT NULL,
	"sex" varchar(1),
	"dob" timestamp NOT NULL,
	"dod" timestamp,
	"hospital_expire_flg" varchar(1) DEFAULT 'N' 
);
CREATE TABLE "d_caregivers" (
	"cgid" integer NOT NULL,
	"label" varchar(6)
);
CREATE TABLE "d_parammap_items" (
	"category" varchar(50) NOT NULL,
	"description" varchar(500)
);
CREATE TABLE "poe_order" (
	"poe_id" bigint NOT NULL,
	"subject_id" integer NOT NULL,
	"hadm_id" integer NOT NULL,
	"icustay_id" integer,
	"start_dt" timestamp,
	"stop_dt" timestamp,
	"enter_dt" timestamp NOT NULL,
	"medication" varchar(255),
	"procedure_type" varchar(50),
	"status" varchar(50),
	"route" varchar(50),
	"frequency" varchar(50),
	"dispense_sched" varchar(255),
	"iv_fluid" varchar(255),
	"iv_rate" varchar(100),
	"infusion_type" varchar(15),
	"sliding_scale" char(1),
	"doses_per_24hrs" float,
	"duration" float,
	"duration_intvl" varchar(15),
	"expiration_val" float,
	"expiration_unit" varchar(50),
	"expiration_dt" timestamp,
	"label_instr" varchar(1000),
	"additional_instr" varchar(1000),
	"md_add_instr" varchar(4000),
	"rnurse_add_instr" varchar(1000)
);
CREATE TABLE "a_meddurations" (
	"subject_id" integer NOT NULL,
	"icustay_id" integer,
	"itemid" integer NOT NULL,
	"elemid" integer NOT NULL,
	"starttime" timestamp NOT NULL,
	"startrealtime" timestamp,
	"endtime" timestamp,
	"cuid" integer,
	"duration" float
);
CREATE TABLE "d_labitems" (
	"itemid" integer NOT NULL,
	"test_name" varchar(50) NOT NULL,
	"fluid" varchar(50) NOT NULL,
	"category" varchar(50) NOT NULL,
	"loinc_code" varchar(7),
	"loinc_description" varchar(100)
);
CREATE TABLE "additives" (
	"subject_id" integer NOT NULL,
	"icustay_id" integer,
	"itemid" integer NOT NULL,
	"ioitemid" integer NOT NULL,
	"charttime" timestamp NOT NULL,
	"elemid" integer NOT NULL,
	"cgid" integer,
	"cuid" integer,
	"amount" float,
	"doseunits" varchar(20),
	"route" varchar(20)
);
CREATE TABLE "poe_med" (
	"poe_id" bigint NOT NULL,
	"drug_type" varchar(20) NOT NULL,
	"drug_name" varchar(100) NOT NULL,
	"drug_name_generic" varchar(100),
	"prod_strength" varchar(255),
	"form_rx" varchar(25),
	"dose_val_rx" varchar(100),
	"dose_unit_rx" varchar(50),
	"form_val_disp" varchar(50),
	"form_unit_disp" varchar(50),
	"dose_val_disp" float,
	"dose_unit_disp" varchar(50),
	"dose_range_override" varchar(2000)
);
CREATE TABLE "totalbalevents" (
	"subject_id" integer NOT NULL,
	"icustay_id" integer,
	"itemid" integer NOT NULL,
	"charttime" timestamp NOT NULL,
	"elemid" integer NOT NULL,
	"realtime" timestamp NOT NULL,
	"cgid" integer,
	"cuid" integer,
	"pervolume" float,
	"cumvolume" float,
	"accumperiod" varchar(20),
	"approx" varchar(20),
	"reset" float,
	"stopped" varchar(20)
);
CREATE TABLE "d_meditems" (
	"itemid" integer NOT NULL,
	"label" varchar(20)
);
CREATE TABLE "db_schema" (
	"created_dt" timestamp DEFAULT LOCALTIMESTAMP ,
	"created_by" varchar(15) DEFAULT USER ,
	"updated_dt" timestamp DEFAULT LOCALTIMESTAMP ,
	"updated_by" varchar(15) DEFAULT USER ,
	"schema_dt" timestamp DEFAULT LOCALTIMESTAMP ,
	"version" varchar(25) NOT NULL,
	"comments" varchar(250)
);
CREATE TABLE "parameter_mapping" (
	"param1_str" varchar(50),
	"param1_num" float,
	"category" varchar(50) NOT NULL,
	"param2_str" varchar(50),
	"param2_num" float,
	"order_num" float,
	"valid_flg" char(1) NOT NULL,
	"comments" varchar(255)
);
CREATE TABLE "a_iodurations" (
	"subject_id" integer NOT NULL,
	"icustay_id" integer,
	"itemid" integer NOT NULL,
	"elemid" integer NOT NULL,
	"starttime" timestamp NOT NULL,
	"startrealtime" timestamp,
	"endtime" timestamp,
	"cuid" integer,
	"duration" float
);
CREATE TABLE "noteevents" (
	"subject_id" integer NOT NULL,
	"hadm_id" integer,
	"icustay_id" integer,
	"elemid" integer,
	"charttime" timestamp NOT NULL,
	"realtime" timestamp,
	"cgid" integer,
	"correction" char(1),
	"cuid" integer,
	"category" varchar(26),
	"title" varchar(255),
	"text" text,
	"exam_name" varchar(100),
	"patient_info" varchar(4000)
);
CREATE TABLE "ioevents" (
	"subject_id" integer NOT NULL,
	"icustay_id" integer,
	"itemid" integer NOT NULL,
	"charttime" timestamp NOT NULL,
	"elemid" integer NOT NULL,
	"altid" integer,
	"realtime" timestamp,
	"cgid" integer,
	"cuid" integer,
	"volume" float,
	"volumeuom" varchar(20),
	"unitshung" float,
	"unitshunguom" varchar(20),
	"newbottle" float,
	"stopped" varchar(20),
	"estimate" varchar(20)
);
CREATE TABLE "labevents" (
	"subject_id" integer NOT NULL,
	"hadm_id" integer,
	"icustay_id" integer,
	"itemid" integer NOT NULL,
	"charttime" timestamp NOT NULL,
	"value" varchar(100),
	"valuenum" float,
	"flag" varchar(10),
	"valueuom" varchar(10)
);
CREATE TABLE "chartevents" (
	"subject_id" integer NOT NULL,
	"icustay_id" integer,
	"itemid" integer NOT NULL,
	"charttime" timestamp NOT NULL,
	"elemid" integer NOT NULL,
	"realtime" timestamp NOT NULL,
	"cgid" integer,
	"cuid" integer,
	"value1" varchar(110),
	"value1num" float,
	"value1uom" varchar(20),
	"value2" varchar(110),
	"value2num" float,
	"value2uom" varchar(20),
	"resultstatus" varchar(20),
	"stopped" varchar(20)
);
CREATE TABLE "procedureevents" (
	"subject_id" integer NOT NULL,
	"hadm_id" integer NOT NULL,
	"itemid" integer,
	"sequence_num" integer NOT NULL,
	"proc_dt" timestamp
);
CREATE TABLE "d_codeditems" (
	"itemid" integer NOT NULL,
	"code" varchar(10),
	"type" varchar(12),
	"category" varchar(13),
	"label" varchar(100),
	"description" varchar(100)
);
CREATE TABLE "d_demographicitems" (
	"itemid" integer NOT NULL,
	"label" varchar(50),
	"category" varchar(19)
);
CREATE TABLE "comorbidity_scores" (
	"subject_id" integer NOT NULL,
	"hadm_id" integer NOT NULL,
	"category" char(10),
	"congestive_heart_failure" float,
	"cardiac_arrhythmias" float,
	"valvular_disease" float,
	"pulmonary_circulation" float,
	"peripheral_vascular" float,
	"hypertension" float,
	"paralysis" float,
	"other_neurological" float,
	"chronic_pulmonary" float,
	"diabetes_uncomplicated" float,
	"diabetes_complicated" float,
	"hypothyroidism" float,
	"renal_failure" float,
	"liver_disease" float,
	"peptic_ulcer" float,
	"aids" float,
	"lymphoma" float,
	"metastatic_cancer" float,
	"solid_tumor" float,
	"rheumatoid_arthritis" float,
	"coagulopathy" float,
	"obesity" float,
	"weight_loss" float,
	"fluid_electrolyte" float,
	"blood_loss_anemia" float,
	"deficiency_anemias" float,
	"alcohol_abuse" float,
	"drug_abuse" float,
	"psychoses" float,
	"depression" float
);
CREATE TABLE "admissions" (
	"hadm_id" integer NOT NULL,
	"subject_id" integer NOT NULL,
	"admit_dt" timestamp NOT NULL,
	"disch_dt" timestamp NOT NULL
);
CREATE TABLE "demographicevents" (
	"subject_id" integer NOT NULL,
	"hadm_id" integer NOT NULL,
	"itemid" integer NOT NULL
);
CREATE TABLE "d_careunits" (
	"cuid" integer NOT NULL,
	"label" varchar(20)
);
CREATE TABLE "icustay_detail" (
	"icustay_id" integer,
	"subject_id" integer,
	"gender" varchar(1),
	"dob" timestamp NOT NULL,
	"dod" timestamp,
	"expire_flg" varchar(1),
	"subject_icustay_total_num" float,
	"subject_icustay_seq" float,
	"hadm_id" integer,
	"hospital_total_num" float,
	"hospital_seq" float,
	"hospital_first_flg" char(1),
	"hospital_last_flg" char(1),
	"hospital_admit_dt" timestamp,
	"hospital_disch_dt" timestamp,
	"hospital_los" float,
	"hospital_expire_flg" char(1),
	"icustay_total_num" float,
	"icustay_seq" float,
	"icustay_first_flg" char(1),
	"icustay_last_flg" char(1),
	"icustay_intime" timestamp NOT NULL,
	"icustay_outtime" timestamp NOT NULL,
	"icustay_admit_age" float,
	"icustay_age_group" varchar(7),
	"icustay_los" float NOT NULL,
	"icustay_expire_flg" char(1),
	"icustay_first_careunit" varchar(20),
	"icustay_last_careunit" varchar(20),
	"icustay_first_service" varchar(110),
	"icustay_last_service" varchar(110),
	"height" float,
	"weight_first" float,
	"weight_min" float,
	"weight_max" float,
	"sapsi_first" float,
	"sapsi_min" float,
	"sapsi_max" float,
	"sofa_first" float,
	"sofa_min" float,
	"sofa_max" float,
	"matched_waveforms_num" float
);
CREATE TABLE "demographic_detail" (
	"subject_id" integer NOT NULL,
	"hadm_id" integer NOT NULL,
	"marital_status_itemid" integer,
	"marital_status_descr" varchar(50),
	"ethnicity_itemid" integer,
	"ethnicity_descr" varchar(50),
	"overall_payor_group_itemid" integer,
	"overall_payor_group_descr" varchar(50),
	"religion_itemid" integer,
	"religion_descr" varchar(50),
	"admission_type_itemid" integer,
	"admission_type_descr" varchar(50),
	"admission_source_itemid" integer,
	"admission_source_descr" varchar(50)
);
CREATE TABLE "a_chartdurations" (
	"subject_id" integer NOT NULL,
	"icustay_id" integer,
	"itemid" integer NOT NULL,
	"elemid" integer NOT NULL,
	"starttime" timestamp NOT NULL,
	"startrealtime" timestamp NOT NULL,
	"endtime" timestamp,
	"cuid" integer,
	"duration" float
);
CREATE TABLE "d_chartitems_detail" (
	"label" varchar(110),
	"label_lower" varchar(110),
	"itemid" integer,
	"category" varchar(50),
	"description" varchar(255),
	"value_type" char(1),
	"value_column" varchar(6),
	"rows_num" float,
	"subjects_num" float,
	"chart_vs_realtime_delay_mean" float,
	"chart_vs_realtime_delay_stddev" float,
	"value1_uom_num" float,
	"value1_uom_has_nulls" char(1),
	"value1_uom_sample1" varchar(20),
	"value1_uom_sample2" varchar(20),
	"value1_distinct_num" float,
	"value1_has_nulls" char(1),
	"value1_sample1" varchar(110),
	"value1_sample2" varchar(110),
	"value1_length_min" float,
	"value1_length_max" float,
	"value1_length_mean" float,
	"value1num_min" float,
	"value1num_max" float,
	"value1num_mean" float,
	"value1num_stddev" float,
	"value2_uom_num" float,
	"value2_uom_has_nulls" char(1),
	"value2_uom_sample1" varchar(20),
	"value2_uom_sample2" varchar(20),
	"value2_distinct_num" float,
	"value2_has_nulls" char(1),
	"value2_sample1" varchar(110),
	"value2_sample2" varchar(110),
	"value2_length_min" float,
	"value2_length_max" float,
	"value2_length_mean" float,
	"value2num_min" float,
	"value2num_max" float,
	"value2num_mean" float,
	"value2num_stddev" float
);
CREATE TABLE "drgevents" (
	"subject_id" integer NOT NULL,
	"hadm_id" integer NOT NULL,
	"itemid" integer NOT NULL,
	"cost_weight" float
);
CREATE TABLE "icustayevents" (
	"icustay_id" integer NOT NULL,
	"subject_id" integer NOT NULL,
	"intime" timestamp NOT NULL,
	"outtime" timestamp NOT NULL,
	"los" float NOT NULL,
	"first_careunit" integer,
	"last_careunit" integer
);
CREATE TABLE "d_chartitems" (
	"itemid" integer NOT NULL,
	"label" varchar(110),
	"category" varchar(50),
	"description" varchar(255)
);
CREATE TABLE "icd9" (
	"subject_id" integer NOT NULL,
	"hadm_id" integer NOT NULL,
	"sequence" integer NOT NULL,
	"code" varchar(100) NOT NULL,
	"description" varchar(255)
);
CREATE TABLE "censusevents" (
	"census_id" integer NOT NULL,
	"subject_id" integer NOT NULL,
	"intime" timestamp NOT NULL,
	"outtime" timestamp NOT NULL,
	"careunit" integer,
	"destcareunit" integer,
	"dischstatus" varchar(20),
	"los" float,
	"icustay_id" integer
);
CREATE TABLE "microbiologyevents" (
	"subject_id" integer,
	"hadm_id" integer,
	"charttime" timestamp,
	"spec_itemid" integer,
	"org_itemid" integer,
	"isolate_num" float,
	"ab_itemid" integer,
	"dilution_amount" varchar(72),
	"dilution_comparison" varchar(10),
	"interpretation" varchar(1)
);
CREATE TABLE "icustay_days" (
	"icustay_id" integer,
	"subject_id" integer,
	"seq" integer,
	"begintime" timestamp,
	"endtime" timestamp,
	"first_day_flg" char(1),
	"last_day_flg" char(1)
);
CREATE TABLE "medevents" (
	"subject_id" integer NOT NULL,
	"icustay_id" integer,
	"itemid" integer NOT NULL,
	"charttime" timestamp NOT NULL,
	"elemid" integer NOT NULL,
	"realtime" timestamp NOT NULL,
	"cgid" integer,
	"cuid" integer,
	"volume" float,
	"dose" float,
	"doseuom" varchar(20),
	"solutionid" integer,
	"solvolume" float,
	"solunits" varchar(20),
	"route" varchar(20),
	"stopped" varchar(20)
);
CREATE TABLE "d_ioitems" (
	"itemid" integer NOT NULL,
	"label" varchar(600),
	"category" varchar(50)
);
CREATE TABLE "deliveries" (
	"subject_id" integer NOT NULL,
	"icustay_id" integer,
	"ioitemid" integer NOT NULL,
	"charttime" timestamp NOT NULL,
	"elemid" integer NOT NULL,
	"cgid" integer,
	"cuid" integer,
	"site" varchar(20),
	"rate" float,
	"rateuom" varchar(20) NOT NULL
);




























































































