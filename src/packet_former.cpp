#include "packer_former.h"
#include "crcinterface.h"

PacketFormer::PacketFormer(
        std::ring<uint8_t> * input_buffer,
        std::ring<Packet> * output_buffer,
        uint32_t sender_address,
        uint32_t reciever_address,
        uint8_t data_length) :
        //inputBuffer(input_buffer),
        //outputBuffer(output_buffer),
        senderAddress(sender_address),
        recieverAddress(reciever_address),
        dataLength(data_length),
        sequenceNumber(0) {

    // Make crc generator with polynomial for ISO 3309 in reversed format
    crcGenny = crcutil_interface::CRC::Create(0xEDB88320, 0, 32, true, 0, 0, 0, true, NULL);

}

void PacketFormer::formNextPacket() {

    static Packet tempPacket;
    crcutil_interface::UINT64 tempCRC;

    // Reset packet contents
    (tempPacket.data).clear();
    
    // Data packet type
    tempPacket.controlCode = 0b0000;

    // Data length
    tempPacket.dataLength = dataLength;

    // Add sender and reciever addresses
    tempPacket.senderAddress = senderAddress;
    tempPacket.recieverAddress = recieverAddress;

    // Add current packet number
    tempPacket.sequenceNumber = sequenceNumber;

    // Form data section
    for(int i = dataLength - 1; i >= 0; i--) {

        (tempPacket.data).push_back(inputBuffer->front());
        inputBuffer->pop();

    }

    // Generate CRC and store to ERC field
    crcGenny->Compute(tempPacket.data.data(), tempPacket.dataLength, &tempCRC, NULL);
    tempPacket.erc = (uint32_t) tempCRC;

    // Push formed packet to output buffer
    outputBuffer->push(tempPacket);

}

void PacketFormer::setDataLength(uint8_t data_length) {

    dataLength = data_length;

}

void PacketFormer::resetSequenceNumber() {

    sequenceNumber = 0;
    
}