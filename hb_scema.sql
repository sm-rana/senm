/* /
set search_path to HB1,public;
DROP TABLE IF EXISTS msmts; 
DROP TABLE IF EXISTS channels; 
-- DROP TABLE IF EXISTS Distros;

CREATE TABLE channels (  
			id		SERIAL, 
			net_id_str	VARCHAR(255) NOT NULL,
			msmt_t	CHAR NOT NULL,
			l_lim	DOUBLE PRECISION,
			r_lim	DOUBLE PRECISION,
			stde	DOUBLE PRECISION,
			info	VARCHAR(255),
			
			PRIMARY KEY (id)
			 
			); 
create index on channels (id);

CREATE TABLE msmts (
			id		SERIAL, 
			time	TIMESTAMP NOT NULL, 
			value	DOUBLE PRECISION NOT NULL, 
			cid		INTEGER NOT NULL,
			PRIMARY KEY (id),
			CONSTRAINT fk_chnl FOREIGN KEY (cid) REFERENCES Channels(id) 
			);

create index on msmts (id);

-- Online sensors - flowrates   Q
INSERT INTO channels (net_id_str, msmt_t, stde) VALUES
	('L11059', 'Q', 50),
	('L5160', 'Q', 10),
	('L7358', 'Q', 10),
	('L7513', 'Q', 15),
	('L3372', 'Q', 110),
	('L3371', 'Q', 70);

-- Online sensors - pressures, P
INSERT INTO channels (net_id_str, msmt_t, stde) VALUES
	('4783', 'P', 3.5),
	('7083', 'P', 4.5),
	('11571', 'P', 4.0),
	('3372', 'P', 3.5),
	('10498', 'P', 3.5),
	('7021', 'P', 3.5),
	('6335', 'P', 4.5);

-- Levels L    
INSERT INTO channels (net_id_str, msmt_t) VALUES
	('11569', 'L'),
	('11570', 'L'),  
	('Riverview-new', 'L'),
	('Bloom-new', 'L') ;  

-- Online sensors - Pump flow rates     F
INSERT INTO channels (net_id_str, msmt_t) VALUES
	('p1322', 'F'), 
	('p13724', 'F'), 
	('p13721', 'F'), 
	('p3117', 'F');

-- Pump pressures   B
INSERT INTO channels (net_id_str, msmt_t) VALUES
	('3349', 'B'),
	('3358', 'B'),
	('3371', 'B');

-- Pump control components   C
INSERT INTO channels (net_id_str, msmt_t) VALUES 
	('p1322', 'C'), 
	('p13724', 'C'), 
	('p13721', 'C'), 
	('p3117', 'C');

-- Valve flow F
INSERT INTO channels (net_id_str, msmt_t) VALUES
	('p1325', 'F'), 
	('p13719', 'F');

-- Valve Pressure B
INSERT INTO channels (net_id_str, msmt_t) VALUES
	('1545', 'B'),
	('3352', 'B');

-- Valve Control C
INSERT INTO channels (net_id_str, msmt_t) VALUES
	('p1325', 'C'), 
	('p13719', 'C');
    
-- Known Demands
INSERT INTO channels (net_id_str, msmt_t) VALUES
	('MillerWTP', 'D'); 


COPY msmts FROM 'D:\Dropbox\Real-Time-Modeling\Paulo\Validation HB\v 1.0 - destiny to crash\hb_data_ready_v1.0.csv' DELIMITER ',';
--D:\Dropbox\Real-Time-Modeling\Rana\Data\v3.2_Paulo\dataReady_v3.2.csv
/ */

SELECT * FROM msmts;
SELECT * FROM channels
