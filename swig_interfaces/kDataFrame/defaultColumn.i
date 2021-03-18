class Column{
        public:
        Column(){}
        virtual ~Column(){}

        virtual Column* getTwin()=0;
        virtual void setSize(uint32_t size)=0;
        static Column* getContainerByName(size_t name);

        virtual void serialize(string filename)=0;
        virtual void deserialize(string filename)=0;

        virtual void setValueFromColumn(Column* Container, uint32_t inputOrder,uint32_t outputOrder){

        }


};